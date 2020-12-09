#!/usr/bin/env python
# coding: utf-8

# # GCB SEEK API demo
# 
# ## Housekeeping and libraries

# This demo is derived from the SEEK API examples on Github: https://github.com/seek4science/seekAPIexamples
# 
# Alan Williams, Hadas Leonov, Stuart Owen

# Import the libraries so that they can be used within the notebook
# 
# * **requests** is used to make HTTP calls
# * **json** is used to encode and decode strings into JSON
# * **string** is used to perform text manipulation and checking
# * **getpass** is used to do non-echoing password input
# * **os** is used later for reading directories

import requests
import json
import os
import argparse


# The **base_url** holds the URL to the SEEK instance that will be used in the notebook
# 
# **headers** holds the HTTP headers that will be sent with every HTTP call
# 
# * **Content-type: application/vnd.api+json** - indicates that any data sent will be in JSON API format
# * **Accept: application/vnd.api+json** - indicates that the notebook expects any data returned to be in JSON API format
# * **Accept-Charset: ISO-8859-1** - indicates that the notebook expects any text returned to be in ISO-8859-1 character set

def main(args):

    base_url = args.url
    #base_url = 'http://193.197.73.91'

    headers = {"Content-type": "application/vnd.api+json",
               "Accept": "application/vnd.api+json",
               "Accept-Charset": "ISO-8859-1"}

    # Create a **requests** HTTP **Session**. A **Session** has re-usable settings such as **headers**
    # 
    # The **authorization** is username and password. The user is prompted for this information.


    session = requests.Session()
    session.headers.update(headers)
    session.auth = (args.user, args.password)


    # This is the project in which all of this is created and done
    # it also important for the visibility of the data

    #Please change to the real id
    containing_project_id = args.project_id #3
    assay_id = args.assay_id #5

    # Upload a file from your disk to SEEK

    # Give the local file you want to upload
    local_data_file_name = args.input

    local_data_file = json.loads("""{
      "data": {
        "type": "data_files",
        "attributes": {
          "title": "%s",
          "description": "%s",
          "tags": [
            "example",
            "GCB"
          ],
          "license": "CC-BY-4.0",
          "other_creators": "Bjoern Gruening",
          "content_blobs": [
            {
              "original_filename": "%s",
              "content_type": "%s"
            }
          ],
          "policy": {
            "access": "download",
            "permissions": [
              {
                "resource": {
                  "id": "2",
                  "type": "projects"
                },
                "access": "download"
              }
            ]
          }
        },
        "relationships": {
          "projects": {
            "data": [
              {
                "id": "%s",
                "type": "projects"
              }
            ]
          },
          "assays": {
            "data": [
              {
                "id": "%s",
                "type": "assays"
              }
            ]
          }
        }
      }
    }""" % (args.title, args.description, args.original_filename, args.content_type, containing_project_id, assay_id))

    # local_data_file
    # This is _almost_ exactly like creating a new remote URL. However, we announce that we will add a blob.
    # 
    # Now we first post that and get the result.


    r = session.post(base_url + '/data_files', json=local_data_file)
    r.raise_for_status()


    populated_local_data_file = r.json()
    local_data_file_url = populated_local_data_file['data']['links']['self']
    # populated_local_data_file
    #print(local_data_file_url)


    # We get back information about the **blob** **URL**. To **that** **URL** we upload the actual blob content. 


    blob_url = populated_local_data_file['data']['attributes']['content_blobs'][0]['link']
    # account for possible SEEK misconfiguration
    blob_url = blob_url.replace("http://localhost:3000", base_url) # replace default localhost url by actual one

    #print("Data file will be uploaded to " + blob_url )


    # open the file
    blob_file_path = local_data_file_name
    blob_file = open(blob_file_path, 'rb')


    # Now, tell the header that the request will be an octet stream.

    file_upload_headers = headers.copy()
    file_upload_headers['Content-type'] = "application/octet-stream"


    # Now put this into the session.
    r = session.put(blob_url,  headers = file_upload_headers, data=blob_file)
    r.raise_for_status()
    full_url = base_url + local_data_file_url
    print("Datafile %s is uploaded to <a href='%s'>%s</a>" % (local_data_file_name, full_url, full_url))

    # Set cleanup=True only if you want to clean up after yourself. 
    # If you want to leave a trace in the database (you normally **do**) then
    # keep it at False.

    cleanup = False # don't clean up
    if cleanup:
        # delete all content added above
        session.delete(base_url + local_data_file_url)
    # Close the HTTP **session**

    session.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Seek uploading tool')
    parser.add_argument('--url', help='Seek url', default="http://193.197.73.91")
    parser.add_argument('--user', help='user name', required=True)
    parser.add_argument('--title', help='Dataset title', default="Example title")
    parser.add_argument('--description', help='Description', default="Example description")
    parser.add_argument('--password', help='password', required=True)
    parser.add_argument('--assay_id', type=int, help='Assay ID', required=True)
    parser.add_argument('--project_id', type=int, help='Project ID', required=True)
    parser.add_argument('--input', help='Input file path', required=True)
    parser.add_argument('--content_type', help='content type', default="text/csv")
    parser.add_argument('--original_filename', help='Original filename', default="data")

    args = parser.parse_args()
    main(args)

