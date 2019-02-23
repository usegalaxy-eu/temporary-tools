import argparse
import numpy as np
import math
import imreg_dft as ird
import pandas as pd
import warnings
from PIL import Image

import skimage.io
import skimage.morphology
import openslide
from skimage import transform as tf
from skimage import filters as fil
from skimage.measure import label, regionprops
from scipy.ndimage import distance_transform_edt
from scipy.ndimage.morphology import binary_closing
from scipy.stats import wasserstein_distance
from scipy.optimize import linear_sum_assignment
from skimage.morphology import convex_hull_image, dilation


def HE_MALDI_registration(imgHE_org, imgMS_org, warp_matrix, tissue_array, binary_img):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        Image.MAX_IMAGE_PIXELS = 50000*50000
        #get HE image
        imgHE_org_os = openslide.open_slide(imgHE_org)
        imgHE_org = np.array(imgHE_org_os.get_thumbnail(imgHE_org_os.dimensions))
        #get MS image
        imgMS_org = skimage.io.imread(imgMS_org)

        #defining functions for histogram and difference calculation
        def angle_vector(p1, p2):
            if ((p2[0]-p1[0])*(p2[1]-p1[1])) != 0 :
                vec_ab = [p2[0] - p1[0], p2[1] - p1[1]]
                vec_norm = [1, 0]
                skalar = np.dot(vec_ab, vec_norm)
                length_ab = math.sqrt(vec_ab[0]*vec_ab[0]+vec_ab[1]*vec_ab[1])
                length_norm = math.sqrt(vec_norm[0]*vec_norm[0] + vec_norm[1]*vec_norm[1])
                argument = skalar/(length_ab*length_norm)
                if argument < -1:
                    argument = -1
                elif argument > 1:
                    argument = 1
                angle = math.acos(argument)
            else:
                angle = 0
            return angle

        def histo_extract(angle_list, bin_size):
            histo = np.histogram(np.asarray(angle_list), bins = bin_size)[0]
            histo_rolled = np.roll(histo, -(np.argmax(histo)))
            return histo_rolled


        def hist_diff_calc (histo1, histo2):
            diff = wasserstein_distance(histo1, histo2)
            return diff

        def distance_calc(p1, p2):
            vec_ab = [p2[0] - p1[0], p2[1] - p1[1]]
            dist = math.sqrt(vec_ab[0]**2+vec_ab[1]**2)
            return dist

        #define function for ICP
        def hungarian_w_outlier_rem(src_pts, dst_pts, bin_size):
            matched_src = np.empty((len(src_pts), 2),dtype= np.float32)
            matched_dst = np.empty((len(src_pts), 2),dtype= np.float32)
            cost_matrix = np.empty((len(src_pts), len(dst_pts)),dtype= np.float32)
            for src_i in range(0,len(src_pts)):
                for dst_i in range(0, len(dst_pts)):
                    cost_matrix[src_i, dst_i] = distance_calc(src_pts[src_i], dst_pts[dst_i])
            matches = linear_sum_assignment(cost_matrix)
            matched_src = src_pts[matches[0]]
            matched_dst = dst_pts[matches[1]]

            histolist_src = np.empty((len(matched_src), 2+bin_size),dtype= np.float32)
            histolist_dst = np.empty((len(matched_dst), 2+bin_size),dtype= np.float32)
            for i in range(0, len(matched_src)):
                histolist_src[i,0:2] = matched_src[i]
                anglelist_src = []
                for j in range(0, len(src_pts)):
                    angle_src = angle_vector(matched_src[i], src_pts[j])
                    anglelist_src.append(angle_src)
                histolist_src[i, 2:] = histo_extract(anglelist_src, bin_size)
            for i in range(0, len(matched_dst)):
                histolist_dst[i,0:2] = matched_dst[i]
                anglelist_dst = []
                for j in range(0, len(dst_pts)):
                    angle_dst = angle_vector(matched_dst[i], dst_pts[j])
                    anglelist_dst.append(angle_dst)
                histolist_dst[i, 2:] = histo_extract(anglelist_dst, bin_size)
            matches_with_histo = np.empty((len(histolist_src), 5), dtype=np.float32)
            for i in range(0, len(histolist_src)):
                matches_with_histo[i,0:2] = histolist_src[i,0:2]
                matches_with_histo[i,2:4] = histolist_dst[i,0:2]
                matches_with_histo[i,4] = hist_diff_calc(histolist_src[i][2:], histolist_dst[i][2:])
            sort_order = np.argsort(matches_with_histo[:,4])
            sorted_matches = matches_with_histo[sort_order]
            #perc_outliers = int(float(sorted_matches.shape[0])*0.4)
            #matches_src_wo_outliers = sorted_matches[0:perc_outliers,0:2]
            #matches_dst_wo_outliers = sorted_matches[0:perc_outliers,2:4]
            threshold_outliers = sorted_matches[int(sorted_matches.shape[0]*0.3),4]
            del_outliers = np.argwhere(list(map(lambda x: x > threshold_outliers, matches_with_histo[:,4])))
            matches_wo_outliers = np.delete(matches_with_histo, del_outliers, axis=0)
            matches_src_wo_outliers = matches_wo_outliers[:,0:2]
            matches_dst_wo_outliers = matches_wo_outliers[:,2:4]
            if matches_src_wo_outliers.shape[0] < 4 or matches_dst_wo_outliers.shape[0] < 4:
                matches_src_wo_outliers = sorted_matches[0:4,0:2]
                matches_dst_wo_outliers = sorted_matches[0:4,2:4]
            return matches_src_wo_outliers, matches_dst_wo_outliers

        if tissue_array == True:
            #openslide downscaling
            long_edge = 5000
            imgHE_org_s_pre = np.array(imgHE_org_os.get_thumbnail((long_edge, long_edge)))
            imgHE_org_s = (imgHE_org_s_pre).astype(np.float32)/255
            factor_openslide = 5000.0/np.max(imgHE_org.shape)
            #remove black lines in HE image
            imgHE_02merged = imgHE_org_s[:,:,0]+(imgHE_org_s[:,:,2])
            RBadd = np.sort(np.ravel(imgHE_02merged))
            percbl = int(RBadd.shape[0]*0.01)
            thresholdbl = RBadd[percbl]
            blrmask = imgHE_02merged < thresholdbl
            imgHE_org_s[blrmask] = np.max(np.reshape(imgHE_org_s, [-1, 3]), axis=0)

            #gaussian filter
            imgHE_gf = skimage.filters.gaussian(imgHE_org_s, sigma = 3, preserve_range = True,)

            #HSV of HE image; threshold is color range
            imgHE_hsv = skimage.color.rgb2hsv(imgHE_gf)
            #threshold determination threshold HE image:
            value_listHE = np.sort(np.ravel(imgHE_hsv[:,:,1]))
            perctrHE = int(value_listHE.shape[0]*0.7)
            thresholdHE = value_listHE[perctrHE]

            if binary_img == False:
                #median filter of green values of MALDI
                imgMS_hsv = skimage.color.rgb2hsv(imgMS_org)
                #med_mask = skimage.morphology.square(9)
                #imgMS_sat = imgMS_hsv[:,:,2]
                imgMS_med = imgMS_hsv[:,:,2]#fil.median(imgMS_sat, med_mask)
                #threshold determination MS
                value_list = np.sort(np.ravel(imgMS_med))
                perctr = int(value_list.shape[0]*0.2)
                thresholdMS = value_list[perctr]
                #segmentation MS
                imgMS_seg = (imgMS_med < thresholdMS).astype(np.float32)*255
            else:
                imgMS_inv = np.invert(imgMS_org)
                imgMS_inv[0:2,:] = 0
                imgMS_inv[-2:,:] = 0
                imgMS_inv[:,0:2] = 0
                imgMS_inv[:,-3:] = 0
                imgMS_seg = (imgMS_inv[:,:,2]>0).astype(np.float32)*255
                    
            #segmentation HE
            imgHE_seg = np.logical_and(imgHE_hsv[:,:,0] > 280.0/360, imgHE_hsv[:,:,1] > thresholdHE).astype(np.float32)*255

            #generate binary images
            imgHE = imgHE_seg/255
            imgMS = imgMS_seg/255
            #downscale HE image to size of MS image
            factor_HE2MS = float(imgMS.shape[1])/imgHE.shape[1]
            imgHE_sz = tf.rescale(imgHE, scale=factor_HE2MS, mode="constant")
            #labeling of HE image
            imgHE_lb_pre_rem = label(imgHE_sz>0)
            #labeling of MS image
            imgMS_lb_pre_rem = label(imgMS>0)

            #avarage area of regions
            area_HE_sort = np.sort(np.array([p.area for p in regionprops(imgHE_lb_pre_rem)]))


            len_areaHE = int(len(area_HE_sort))
            lenperc_areaHE = int(len(area_HE_sort)*0.5)
            #area_HE_avg = np.mean(area_HE_sort[lenperc_areaHE:len_areaHE])/2
            area_HE_avg = np.max(area_HE_sort)/3
            #removal of small regions HE
            imgHE_lb = imgHE_sz.copy()
            for prop in regionprops(imgHE_lb_pre_rem):
                if prop.area < area_HE_avg:
                    imgHE_lb[imgHE_lb_pre_rem==prop.label] = 0

            area_MS_sort = np.sort(np.array([p.area for p in regionprops(imgMS_lb_pre_rem)]))
            area_MS_avg = np.max(area_MS_sort)/3
            imgMS_lb = imgMS.copy()
            for prop in regionprops(imgMS_lb_pre_rem):
                if prop.area < area_MS_avg:
                    imgMS_lb[imgMS_lb_pre_rem==prop.label] = 0
            #second labeling with removed areas and amount of features MS
            imgHE_lb2 = label(imgHE_lb > 0)
            imgMS_lb2 = label(imgMS_lb > 0)
            con_NrHE = len(regionprops(imgHE_lb2))
            con_NrMS = len(regionprops(imgMS_lb2))

            #extraction of features (centroids)
            propsHE = np.empty((con_NrHE, 2), dtype= np.float32)
            propsMS = np.empty((con_NrMS, 2), dtype= np.float32)
            HE_counter = 0
            MS_counter = 0
            for prop in regionprops(imgHE_lb2):
                propsHE[HE_counter, 0] = prop.centroid[0]
                propsHE[HE_counter, 1] = prop.centroid[1]
                HE_counter = HE_counter + 1
            for prop in regionprops(imgMS_lb2):
                propsMS[MS_counter, 0] = prop.centroid[0]
                propsMS[MS_counter, 1] = prop.centroid[1]
                MS_counter = MS_counter + 1
            #generate point sets with different rotations
            centroidMS_base_r =regionprops(imgMS_lb.astype(np.int_))[0].centroid
            centroidMS_base = np.concatenate([centroidMS_base_r[1:2],centroidMS_base_r[0:1]],axis=0)
            propsMS_r = np.concatenate([propsMS[:,1:2],propsMS[:,0:1]],axis=1)
            radians = [0, math.pi/2, math.pi, (3*math.pi)/2] #radians for different rotations
            rotations = np.empty((len(radians), 3, 3), dtype=np.float32)
            rotated_coordMS = np.empty((len(radians), len(propsMS_r) ,2), dtype=np.float32)
            for i in range(0, len(rotations)):
                rot = tf.SimilarityTransform(rotation = radians[i]).params
                centroidMS_rot = tf.matrix_transform(centroidMS_base, rot)[0]
                cen_diff_x = centroidMS_base[0]-centroidMS_rot[0] 
                cen_diff_y = centroidMS_base[1]-centroidMS_rot[1]
                rot_and_trans = tf.SimilarityTransform(rotation=radians[i],translation = (cen_diff_x,cen_diff_y)).params
                rotations[i] = rot_and_trans
                rotated_coordMS_r = tf.matrix_transform(propsMS_r, rot_and_trans)
                rotated_coordMS[i] = np.concatenate([rotated_coordMS_r[:,1:2],rotated_coordMS_r[:,0:1]],axis=1)
            match_ind = np.empty((len(propsMS),2), dtype=np.float32)
            counter = 0
            for i in range(0,len(propsMS)):
                match_ind[i,0] = counter
                match_ind[i,1] = counter
                counter += 1
            #ICP for every rotated point set
            matrices_icp = np.empty((len(radians), 3, 3), dtype=np.float32)
            for i in range(0, len(rotations)):
                maxiterations = 15
                NN_coordsHE = propsHE
                NN_coordsMS = rotated_coordMS[i]
                previous_matrix_NN = tf.SimilarityTransform().params
                for itr in range (0, maxiterations):
                    matched_coordsMS, matched_coordsHE = hungarian_w_outlier_rem(NN_coordsMS, NN_coordsHE, 8)
                    matched_coordsHE = np.concatenate([matched_coordsHE[:,1:2],matched_coordsHE[:,0:1]],axis=1)
                    matched_coordsMS = np.concatenate([matched_coordsMS[:,1:2],matched_coordsMS[:,0:1]],axis=1)
                    matrix_NN = tf.estimate_transform("Similarity", matched_coordsMS, matched_coordsHE)
                    new_matrix_NN = matrix_NN.params
                    previous_matrix_NN = np.matmul(new_matrix_NN, previous_matrix_NN)        
                    propsMS_p = np.concatenate([rotated_coordMS[i][:,1:2],rotated_coordMS[i][:,0:1]],axis=1)
                    NN_coordsMS = tf.matrix_transform(propsMS_p, previous_matrix_NN)
                    NN_coordsMS = np.concatenate([NN_coordsMS[:,1:2],NN_coordsMS[:,0:1]],axis=1)
                matrices_icp[i] = np.matmul(previous_matrix_NN, rotations[i])
            #select best registration based on area overlap of template and transformed image
            area_list = np.empty((len(matrices_icp), 1), dtype=np.float32)
            for i in range(0,len(rotations)):
                warped_img = tf.warp(imgMS_lb, np.linalg.inv(matrices_icp[i]), output_shape=(imgHE_lb.shape[0],imgHE_lb.shape[1]))
                diff_img = np.sqrt((imgHE_lb-warped_img)**2)
                area = regionprops(diff_img.astype(np.int_))[0].area
                area_list[i] = area
            pos_best_reg = np.argmin(area_list)
            icp_reg_matrix = matrices_icp[pos_best_reg]
            #fuse all matrices
            scale1 = tf.SimilarityTransform(scale=1/factor_openslide).params
            scale2 = tf.SimilarityTransform(scale=1/factor_HE2MS).params
            scale1_2 = np.matmul(scale2, scale1)
            full_matrix = np.matmul(scale1_2, icp_reg_matrix)
        
        else:
            #remove black lines in HE image
            imgHE_02merged = imgHE_org[:,:,0]+(imgHE_org[:,:,2])
            RBadd = np.sort(np.ravel(imgHE_02merged))
            percbl = int(RBadd.shape[0]*0.01)
            thresholdbl = RBadd[percbl]
            blrmask = imgHE_02merged < thresholdbl
            imgHE_org[blrmask] = np.max(np.reshape(imgHE_org, [-1, 3]), axis=0)

            #gaussian filter
            imgHE_gf = skimage.filters.gaussian(imgHE_org, sigma = 3, preserve_range = True,)

            #threshold determination HE
            imgHE_hsv = skimage.color.rgb2hsv(imgHE_gf)

            #threshold determination threshold HE image:
            value_listHE = np.sort(np.ravel(imgHE_hsv[:,:,1]))
            perctrHE = int(value_listHE.shape[0]*0.7)
            thresholdHE = value_listHE[perctrHE]

            #threshold determination MS
            imgMS_hsv = skimage.color.rgb2hsv(imgMS_org)
            value_list = np.sort(np.ravel(imgMS_hsv[:,:,2]))
            perctr = int(value_list.shape[0]*0.18)
            thresholdMS = value_list[perctr]

            #segmentation HE
            imgHE_seg = np.bitwise_and(imgHE_hsv[:,:,0] > 280/360, imgHE_hsv[:,:,1] > thresholdHE).astype(np.float32)*255
            #segmentation MS
            imgMS_seg = (imgMS_hsv[:,:,2] < thresholdMS).astype(np.float32)*255

            #generate binary images
            imgHE = imgHE_seg/255
            imgMS_re = imgMS_seg/255
            #closing of object in MS image
            closingMS1 = binary_closing(imgMS_re, iterations=1)
            imgMS = imgMS_re
            #labeling of HE image
            imgHE_lb_pre_rem = label(imgHE>0)
            #labeling of MS image
            imgMS_lb_pre_rem = label(imgMS>0)

            #avarage area of regions
            area_HE_sort = np.sort(np.array([p.area for p in regionprops(imgHE_lb_pre_rem)]))
            area_MS_sort = np.sort(np.array([p.area for p in regionprops(imgMS_lb_pre_rem)]))
            area_HE_max = area_HE_sort.max()
            area_MS_max = area_MS_sort.max()

            #selection of object HE
            imgHE_lb = imgHE.copy()
            for prop in regionprops(imgHE_lb_pre_rem):
                if prop.area < area_HE_max:
                    imgHE_lb[imgHE_lb_pre_rem==prop.label] = 0
            #selection of object MS
            imgMS_lb = imgMS.copy()
            for prop in regionprops(imgMS_lb_pre_rem):
                if prop.area < area_MS_max:
                    imgMS_lb[imgMS_lb_pre_rem==prop.label] = 0 
            #generating convex hull image
            convHullHE = convex_hull_image(imgHE_lb)
            convHullMS = convex_hull_image(imgMS_lb)
            #estimating zoom for registration
            axisHE = regionprops(convHullHE.astype(np.int_))[0].major_axis_length 
            axisMS = regionprops(convHullMS.astype(np.int_))[0].major_axis_length 
            sz = convHullMS.shape
            center_imgMS = [int(sz[0]/2), int(sz[1]/2)]
            center_regionHE = regionprops(convHullHE.astype(np.int_))[0].centroid
            scalefactor_HE2MS = axisMS/axisHE
            mat_temp = tf.SimilarityTransform(scale=scalefactor_HE2MS).params
            temp_imgHE = tf.warp(convHullHE, np.linalg.inv(mat_temp), output_shape=(sz[0],sz[1]))
            temp_centerHE = regionprops(temp_imgHE.astype(np.int_))[0].centroid
            diff_y = center_imgMS[0] - int(temp_centerHE[0])
            diff_x = center_imgMS[1] - int(temp_centerHE[1])
            scale_m = tf.SimilarityTransform(translation=(diff_x,diff_y), scale=scalefactor_HE2MS).params
            imgHE_lb_sc = tf.warp(convHullHE, np.linalg.inv(scale_m), output_shape=(sz[0],sz[1]))
            imgHE_lb_sc_dil = dilation(imgHE_lb_sc)
            imgMS_dist = distance_transform_edt(convHullMS)
            imgHE_dist = distance_transform_edt(imgHE_lb_sc_dil)
            #registration with HE as template
            const = {"angle":(0,70)}
            values_reg = ird.imreg.similarity(imgHE_dist, imgMS_dist, numiter=20, order=3, constraints= const, reports=None)
            #generation of registrated matrix
            reg_scale = values_reg["scale"]
            reg_rot = -math.radians(values_reg["angle"])
            reg_trans = values_reg["tvec"]            
            #generation of registrated matrix
            reg_scale = values_reg["scale"]
            reg_rot = -math.radians(values_reg["angle"])
            reg_trans = values_reg["tvec"]

            reg_matrix_s0 = tf.SimilarityTransform(translation=(-(imgMS_dist.shape[1]+1)/2, -(imgMS_dist.shape[0]+1)/2)).params
            reg_matrix_s1 = tf.SimilarityTransform(scale=reg_scale).params
            reg_matrix_s2 = tf.SimilarityTransform(translation=((imgMS_dist.shape[1]+1)/2, (imgMS_dist.shape[0]+1)/2)).params
            reg_matrix_s = np.matmul(reg_matrix_s2, np.matmul(reg_matrix_s1, reg_matrix_s0))

            reg_matrix_0 = tf.SimilarityTransform(translation=(-np.round(imgMS_dist.shape[1]/2), -np.round(imgMS_dist.shape[0]/2))).params
            reg_matrix_1 = tf.SimilarityTransform(rotation=reg_rot).params
            reg_matrix_2 = tf.SimilarityTransform(translation=(np.round(imgMS_dist.shape[1]/2)-0.5,np.round(imgMS_dist.shape[0]/2)-0.5)).params
            reg_matrix_r = np.matmul(reg_matrix_2, np.matmul(reg_matrix_1, reg_matrix_0))

            reg_matrix_t = tf.SimilarityTransform(translation=(reg_trans[1],reg_trans[0])).params
            
            reg_matrix = np.matmul(reg_matrix_t, np.matmul(reg_matrix_r, reg_matrix_s))
            #generation of complete transformation matrix
            sz_HE = imgHE_org.shape
            full_matrix =  np.matmul(np.linalg.inv(scale_m), reg_matrix)

        # convert maldi->he to he->maldi
        full_matrix = np.linalg.inv(full_matrix)

        #save warp matrix as csv file
        pd.DataFrame(full_matrix).to_csv(warp_matrix, header=None, index=False, sep="\t")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Registrate two images")
    parser.add_argument("imgHE_org", help="HE image")
    parser.add_argument("imgMS_org", help="MALDI image")
    parser.add_argument("warp_matrix", help="Output: 3x3 csv matrix")
    parser.add_argument('--tissue_array', dest='tissue_array', action='store_true')
    parser.add_argument('--binary_img', dest='binary_img', action='store_true')
    args = parser.parse_args()
    HE_MALDI_registration(args.imgHE_org, args.imgMS_org, args.warp_matrix, args.tissue_array, args.binary_img)
