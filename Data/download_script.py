#!/usr/bin/env python
  
# -*- coding: utf-8 -*-
#
# Script to loop on timespan (day / week / month or year) to optimize dataset requests (heavy in terms of number of files to be manipulated).


# this script should be run from inside the Data directory
  
import numpy as np
import pandas as pd
import os
import platform
import subprocess
import datetime as dt
import time
import calendar
import glob
import argparse
# Debug
import ipdb
  

def main(year):
    # These have to be passed as inputs
    # product = '030'
    # year = 2008

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # General Parameters - Tools - Proxy Network - Output Directory
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      
    # Module declaration to the motu-client opensource-TOOLS to connect to MOTU CopernicusMarineHub.
    # If you can't call it as module, then input the 'motu-client.py' absolute path. By default, usually in "Downloads" dir. after having followed the article on "python basic requirements":
    # http://marine.copernicus.eu/faq/what-are-the-motu-and-python-requirements/?idpage=169
    # Deprecated : motu_cl = '{absolute_path_to}/motu-client-python/motu-client.py'
    # Deprecated : motu_cl = 'python -m motu-client'
    motu_cl = 'python -m motuclient'
      
    # Copernicus Marine API Key - Login Credentials 
    # To create an account reach: http://marine.copernicus.eu/services-portfolio/register-now/.
    # If already created but forgotten reach: http://marine.copernicus.eu/faq/forgotten-password/?idpage=169
    username_cmems = 'fnencioli'
    password_cmems = 'Francesco@CMEMS2019'
      
    #------------------------------------
    # Proxy Configuration
    # Please replace "False" by "True" if you use a proxy to connect to internet and fill in the below variables.
    proxy_flag = False
    if proxy_flag:
        proxy_server_url = "http://your_proxy_url.com"
        proxy_server_port = "port"
        proxy_user_login = "your_proxy_user_login"
        proxy_user_password = "your_proxy_user_password"
    #------------------------------------
      
    # Output directory name to store the Copernicus Marine data - (do not use whitespace character)
    # If only 'folder-name' is given (not in absolute path), then it will be converted automatically into '$HOME/folder-name/'
    local_storage_directory_name = './'
      
    # - - - - - - - - - - - - - - - - - - - - - - - - -
    # Product(s), Dataset(s) and MOTU server Parameters 
    # - - - - - - - - - - - - - - - - - - - - - - - - -
      
    # CMEMS Variables & Dataset ID & Service ID & MOTU server ID
    # Define a dict to get required parameters of our daily temperature data request. 
    # It should looks like:
    #                      {file_name (defined by yourself for practical reason): \
    #                        [variable (aka -v), \
    #                        dataset_id (aka -d), \
    #                        product_id (aka -s), \
    #                        motu_id (aka -m)]
    #                        }
      
    #  -v VARIABLE
    #  --variable=VARIABLE
    #                        The variable name or standard_name (list of strings, e.g. --variable=thetao or -v sea_water_potential_temperature)
    #  -d PRODUCT_ID
    #  --product-id=PRODUCT_ID
    #                        The product (data set) to download (string e.g. -d global-analysis-forecast-phy-001-024)
    #  -s SERVICE_ID
    #  --service-id=SERVICE_ID
    #                        The service identifier (string e.g. --service=GLOBAL_ANALYSIS_FORECAST_PHY_001_024-TDS)
    #  -m MOTU
    #  --motu=MOTU
    #                        The motu server to use (url, e.g. -m http://nrt.cmems-du.eu/motu-web/Motu or --motu http://my.cmems-du.eu/motu-web/Motu)
      
    # /!\ all CMEMS products are NOT hosted by a single server - they are grouped by MultiYear or NearRealTime products respectively on http://my.cmems-du.eu/motu-web/Motu and http://nrt.cmems-du.eu/motu-web/Motu
    # You can always rely on the "VIEW SCRIPT" button of the Copernicus Marine Website (marine.copernicus.eu),
    # using its DataExtraction WebInterface (also called GUI which stands for Graphical User Interface).
    # It will generate the parameters of your extraction settings based on your selection.
    # Please refer to this article to understand how to call/trigger this webservice/feature: http://marine.copernicus.eu/faq/how-to-write-and-run-the-script-to-download-cmems-products-through-subset-or-direct-download-mechanisms/?idpage=169
     
    product_id = 'dataset-duacs-rep-global-merged-allsat-phy-l4'
    service_id = 'SEALEVEL_GLO_PHY_L4_REP_OBSERVATIONS_008_047-TDS'
    ficout_tags = ['JFM','AMJ','JAS','OND']
    variables = "--variable ugos --variable vgos"
       

 
    dict_id = {"Currents MERCATOR model": \
               ["%s" % variables, \
                "-d %s" % product_id, \
                "-s %s" % service_id, \
                "-m  http://my.cmems-du.eu/motu-web/Motu"]
              }
      
    # - - - - - - - - - - - - - - - - - - - - - -
    # Geographical Area Parameters and Timerange
    # - - - - - - - - - - - - - - - - - - - - - -
      
    #  -y LATITUDE_MIN
    #  --latitude-min=LATITUDE_MIN
    #                        The min latitude (float in the interval [-90 ; 90])
    #  -Y LATITUDE_MAX
    #  --latitude-max=LATITUDE_MAX
    #                        The max latitude (float in the interval [-90 ; 90])
    #  -x LONGITUDE_MIN
    #  --longitude-min=LONGITUDE_MIN
    #                        The min longitude (float in the interval [-180 ; 180])
    #  -X LONGITUDE_MAX
    #  --longitude-max=LONGITUDE_MAX
    #                        The max longitude (float in the interval [-180 ; 180])
    #  -z DEPTH_MIN
    #  --depth-min=DEPTH_MIN
    #                        The min depth (float in the interval [0 ; 2e31] or
    #                        string 'Surface')
    #  -Z DEPTH_MAX
    #  --depth-max=DEPTH_MAX
    #                        The max depth (float in the interval [0 ; 2e31] or
    #                        string 'Surface')
    #  -t DATE_MIN
    #  --date-min=DATE_MIN
    #                        The min date with mandatory hour resolution (string following format "YYYY-MM-DD HH:MM:SS"),
    #                        e.g. -t "2016-06-10 12:00:00". Be careful to NOT forget double quotes around the date.
    #  -T DATE_MAX
    #  --date-max=DATE_MAX
    #                        The max date with mandatory hour resolution (string following format "YYYY-MM-DD HH:MM:SS"),
    #                        e.g. -T "2016-06-10 12:00:00". Be careful to NOT forget double quotes around the date.
      
     
    # Area : x east-west longitude, y north-south latitude, z depth
     
    xmin_longitude = "-180"
    xmax_longitude = "180"
    ymin_latitude = "-90"
    ymax_latitude = "90"
    #if product == '025':
    #    zmin_depth = "13.99"
    #    zmax_depth = "13.9912"
    #elif product == '030':
    #    zmin_depth = "13.467"
    #    zmax_depth = "13.4673"
      
    # Date - Timerange
    # Nencio: Process one year at the time
    yyyystart = year
    yyyyend = year
    mmstart = 1
    mmend = 12
    #hhstart = " 12:00:00"
    hhstart = " 00:00:00" # this was changed to ensure the first day of the 15-day interval is also loaded. I am not sure if this different hour (00:00:00 instead of 12:00:00) will cause trouble later
    hhend = " 12:00:00"
    # Use pandas to create time intervals  
#    if product == '025':
#        dend = pd.date_range(dt.datetime(yyyystart,mmstart,1),dt.datetime(yyyyend,mmend,31),freq='Q')
#        dini = pd.date_range(dt.datetime(yyyystart,mmstart,1),dt.datetime(yyyyend,mmend,31),freq='QS')
#    elif product == '030':
#        dini = pd.date_range(dt.datetime(yyyystart,mmstart,1),dt.datetime(yyyyend,mmend,31),freq='SMS')
#        # Extract values (since datetimeindex cannot be modified)
#        _dini = dini.values
#        # Add one day to 15 of the month
#        _dini[1::2] += np.timedelta64(1,'D')
#        # Convert back to datetimeindex
#        dini = pd.to_datetime(_dini)
#        dend = pd.date_range(dt.datetime(yyyystart,mmstart,1),dt.datetime(yyyyend,mmend,31),freq='SM')
    dini = pd.date_range(dt.datetime(yyyystart,mmstart,1),dt.datetime(yyyyend,mmend,31),freq='SMS')
    # Extract values (since datetimeindex cannot be modified)
    _dini = dini.values
    # Add one day to 15 of the month
    _dini[1::2] += np.timedelta64(1,'D')
    # Convert back to datetimeindex
    dini = pd.to_datetime(_dini)
    dend = pd.date_range(dt.datetime(yyyystart,mmstart,1),dt.datetime(yyyyend,mmend,31),freq='SM')
    
    # Output prefix file name
    # pre_name = "CMEMS_"

      
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #                     Main Program
    #
    #          Motu Client Call through Python Loop
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Specific comment For WINDOWS USER:
    # If you're using this script for the first time, you
    # shouldn't be worried by lines below. In your text editor,
    # just save your script (ctrl + s), quit (alt + F4) and launch it
    # (WinKey + R then input "cmd" then click on ENTER)
    # to finally type in your terminal
    # "python.exe script_name.py"
    #
    # For users, be careful if you have to modify the lines below.
    # CMEMS Central Service Desk will be happy to help you
    # either via email (servicedesk.cmems@mercator-ocean.eu)
    # or via the CMEMS Forum (http://bit.ly/1L1Iy5f)
      
    #-------------------------
    # Nencio: not needed
    #         Use local_storage_directory_name directly
    #
    # # Check if output directory is well formated and if it exists, otherwise create it
    # absolute_path_substring = ['/home/', 'C:\\']
    # if local_storage_directory_name[-1] != '/':
    #     local_storage_directory_name = local_storage_directory_name + "/"
    # if not any(x in local_storage_directory_name for x in absolute_path_substring):
    #     local_storage_directory_name = os.path.expanduser('~') + "/" + local_storage_directory_name
    # if not os.path.exists(local_storage_directory_name):
    #     os.makedirs(local_storage_directory_name)
    #-------------------------
      
    # Flags to let the server clears the buffer - better to be respectful when retrieving OPEN data
    buffer_flag = False
    cmd_flag = False
      
    # Error Handle on dates (to illustrate an if statement with eval param '>')
    if yyyystart > yyyyend:
        print("[ERROR] in [Date Parameters]")
        print("""Please double check your date parameters, specifically the "yyyystart" which is currently greater than "yyyyend.""")
        print("""End of data extraction service.""")
        raise SystemExit
      
    # Other variable definitions to be compatible with deprecated script versions still available on the Internet
    log_cmems = "-u " + username_cmems
    pwd_cmems = "-p " + password_cmems
    # Moved to within the for loop
    # pre_fic_cmd = "-f "+ pre_name
    out_cmd = "-o " + local_storage_directory_name
    #--------------------------
    if proxy_flag:
        proxy_user = "--proxy-user " + proxy_user_login
        proxy_pwd = "--proxy-pwd " + proxy_user_password
        proxy_server = "--proxy-server " + proxy_server_url + ":" + proxy_server_port
    #--------------------------
    xmin = "-x " + xmin_longitude
    xmax = "-X " + xmax_longitude
    ymin = "-y " + ymin_latitude
    ymax = "-Y " + ymax_latitude
    #zmin = "-z " + zmin_depth
    #zmax = "-Z " + zmax_depth
      
    # # To illustrate a simple Error Handle to delete a file when desired
    # try:
    #     os.remove(out_cmd.split()[1] + logfile)
    # except OSError:
    #     print("")

    # Check number of files already downloaded
    flist = glob.glob(os.path.join(local_storage_directory_name,'%s_*_DATASET-DUACS-REP-GLOBAL-MERGED-ALLSAT-PHY-L4.nc' % (year)))
    if len(flist) == len(dini):
        print("All files already downloaded")
        return
    elif len(flist) > len(dini):
        print("Something wrong more files downloaded than available")
        return

     
    print("\n+----------------------------+\n| ! - CONNEXION TO CMEMS HUB |\n+----------------------------+\n\n")
          
    # To illustrate a For_Loop in order to generate download requests for several datasets held in a product
    while len(flist)<len(dini):
        for key, value in dict_id.items():

            # File to log unsuccessful data extraction request(s)
            logfile = 'logfile_%s.txt' % dt.datetime.now().strftime('%Y%m%dT%H%M')

            for i, idini in enumerate(dini):
#                if product == '025':
#                    ficout = '%04d_%s_GLOBAL_REANALYSIS_PHY_001_%s_uv_global.nc' % (year,ficout_tags[i],product)
#                elif product == '030':
#                    ficout = '%04d_%02d%1d_GLOBAL_REANALYSIS_PHY_001_%s_uv_global.nc' % (year,idini.month,2-idini.day%16,product)
                ficout = '%04d_%02d%1d_DATASET-DUACS-REP-GLOBAL-MERGED-ALLSAT-PHY-L4.nc' % (year,idini.month,2-idini.day%16)
               
                if buffer_flag:
                    # print("Little pause to let the server clearing the buffer, it will AUTOMATICALLY resume once it's completed.\nNot mandatory but server-friendly <span class="Emoticon Emoticon1"><span>:-)</span></span>\n")
                    print("Little pause to let the server clearing the buffer, it will AUTOMATICALLY resume once it's completed.")
                    time.sleep(2)
                    buffer_flag = False
                          
                # Date declaration
                date_start = dini[i]
                date_end = dend[i]
                  
                # date_end_cmd = (dt.datetime(date_start.year, date_start.month,\
                #     calendar.monthrange(date_start.year, date_start.month)[1]))
                date_cmd = ' -t \"' + date_start.strftime("%Y-%m-%d") + hhstart + '\"'\
                +' -T \"' + date_end.strftime("%Y-%m-%d") + hhend + '\"'
                # fic_cmd = pre_fic_cmd + key + "_" + date_end_cmd.strftime("%Y-%m") + ".nc"
                fic_cmd = '-f %s' % ficout 
                # ficout = pre_name + key + "_" + date_end_cmd.strftime("%Y-%m") + ".nc"
                print("----------------------------------\n- ! - Processing dataset request : %s"%ficout)
                print("----------------------------------\n")
                # Download only if not already downloaded
                if not os.path.exists(out_cmd.split()[1] + ficout):
                    if proxy_flag:
                        if 'zmin_depth' not in locals():
                            cmd = ' '.join([motu_cl, log_cmems, pwd_cmems,\
                                        value[3], value[2], value[1],\
                                        xmin, xmax, ymin, ymax,\
                                        date_cmd, value[0], out_cmd, fic_cmd,\
                                        proxy_server, proxy_user, proxy_pwd, "-q"])
                        else:
                            cmd = ' '.join([motu_cl, log_cmems, pwd_cmems,\
                                        value[3], value[2], value[1],\
                                        xmin, xmax, ymin, ymax, zmin, zmax,\
                                        date_cmd, value[0], out_cmd, fic_cmd,\
                                        proxy_server, proxy_user, proxy_pwd, "-q"])
                    else:
                        if 'zmin_depth' not in locals():
                            cmd = ' '.join([motu_cl, log_cmems, pwd_cmems,\
                                        value[3], value[2], value[1],\
                                        xmin, xmax, ymin, ymax,\
                                        date_cmd, value[0], out_cmd, fic_cmd, "-q"])
                        else:
                            cmd = ' '.join([motu_cl, log_cmems, pwd_cmems,\
                                        value[3], value[2], value[1],\
                                        xmin, xmax, ymin, ymax, zmin, zmax,\
                                        date_cmd, value[0], out_cmd, fic_cmd, "-q"])
                    print("## MOTU API COMMAND ##")
                    print(cmd)
                    print("\n[INFO] CMEMS server is checking both your credentials and command syntax. If successful, it will extract the data and create your dataset on the fly. Please wait. \n")
                    subpro=subprocess.Popen(cmd,shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                    message,erreur = subpro.communicate()
                    stat = subpro.returncode
                    if stat != 0:
                            print("-- ERROR Incorrect Credentials :\n %s"%message)
                            with open(out_cmd.split()[1] + logfile,'a') as mylog:
                                mylog.write("Error : %s NOK\nDue to : %s"%(ficout,message))
                            print("""[INFO] Failed data extraction has been logged.\n""")
                            if b'HTTP Error 400' in message:
                                print('''[INFO] Copernicus Marine USERNAME ('username_cmems') and/or PASSWORD ('password_cmems') are incorrect.\n\n[INFO] To execute the MOTU API COMMAND from your shell/terminal, please note the following rules:\n
                                On *nix OS, you must use the single quote, otherwise it may expand special characters.
                                [...] -u 'string' or --user='string' [...]\n
                                On Windows OS, you must use the double quote, because single quotes are treated literally.
                                [...] -p "string" or --pwd="string" [...]\n''')
                                raise SystemExit
                            if b'HTTP Error 407' in message:
                                print('''[INFO] Proxy Authentication Required to connect to the Central Authentication System https://cmems-cas.cls.fr/cas/login\n\n[INFO] Check the value of proxy_flag (it should be True).\n\n[INFO] Double check your proxy settings:\n  --proxy-server=PROXY_SERVER\n                        the proxy server (url)\n  --proxy-user=PROXY_USER\n                        the proxy user (string)\n  --proxy-pwd=PROXY_PWD\n                        the proxy password (string)\n\n[INFO] If your proxy credentials are correct but your proxy password (string) contains a '@' then replace it by '%%40' ''')
                                print('''[INFO] This issue is raised due either a misconfiguration in proxy settings or a network issue. If it persists, please contact your network administrator.''')
                                raise SystemExit
                            if b'HTTP Error 403' in message:
                                print('''[INFO] Copernicus Marine USERNAME ('username_cmems') has been suspended.\n[INFO] Please contact our Support Team either:\n  - By mail: servicedesk.cmems@mercator-ocean.eu or \n  - By using a webform, reaching the marine.copernicus.eu website and triggering the ANY QUESTIONS? button.''')
                                raise SystemExit
                    else:
                        if b"[ERROR]" in message:
                            print("-- ERROR Downloading command :\n %s"%message)
                            with open(out_cmd.split()[1] + logfile,'a') as mylog:
                                mylog.write("Error : %s NOK\nDue to : %s"%(ficout,message))
                            print("""[INFO] Failed data extraction has been logged.\n""")
                        else:
                                print("-- MOTU Download successful :\n %s OK\n"%fic_cmd.split()[1])
                                cmd_flag = True
                else:
                    print("-- Your dataset for %s has already been downloaded in %s --\n"% (fic_cmd.split()[1],out_cmd.split()[1]))
                    cmd_flag = False
                  
                if cmd_flag:
                    buffer_flag = True
                    cmd_flag = False
     
        if not os.path.exists(out_cmd.split()[1]+logfile):
            print("\n------------------------------------------------\n - ! - Your Copernicus Dataset(s) are located in %s\n------------------------------------------------\n"%(out_cmd.split()[1]))
        else :
            print("## [ERROR] ##")
            print ("/!\\ Some download requests failed. Please see recommendation in %s%s"%(out_cmd.split()[1], logfile))
        print("+--------------------------------------------+\n| ! - CONNEXION TO CMEMS HUB HAS BEEN CLOSED |\n+--------------------------------------------+\n")
        flist = glob.glob(os.path.join(local_storage_directory_name,'%s_*_DATASET-DUACS-REP-GLOBAL-MERGED-ALLSAT-PHY-L4.nc' % (year)))
    print("\nAll %04d files have been succesfully downloaded" % year)


if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--year', required=True, help="year to process")
    args = parser.parse_args()
    year = int(args.year)
    main(year)
