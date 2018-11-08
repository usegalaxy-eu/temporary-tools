from cairosvg import svg2png
import argparse


def convert(image_svg, image_png):
    svg2png(file_obj = open(image_svg, "rb"), write_to = image_png)
    return image_png

if __name__=="__main__":
    parser = argparse.ArgumentParser(description = "Convert image.svg to image.png")
    parser.add_argument("image_svg", help = "Paste path to image.svg")
    parser.add_argument("image_png", help = "Paste path to directory in which converted image.png should be saved")
    args = parser.parse_args()    
    convert(args.image_svg, args.image_png)

