#!/usr/bin/python

# This function will call into the MODAPS web services server
# and retrieve a space delimited list of URLs to retrieve, which
# will be returned to the invoking shell if called using the bash
# syntax urls=$(python automodis.py)

from SOAPpy import SOAPProxy
import os
import sys
import time
import pdb

def parse_args(args_in):
    """
    Parses command line arguments given in bash. Assumes that all arguments are flag-value pairs
     using long-option nomenclature (e.g. --products MYD06_L2).
    Allows for two types of arguments: required and optional.  Required are specified in the list
     below (req_args) and will cause an error if not present. Optional arguments are given in a
     dictionary; these are arguments that if not specified have reasonable default values, at least
     for implementation for BEHR. The key will be the flag name (without the --) and the value is
     the default value.
    :param args_in: Pass
    :return:
    """
    req_args = ["products", "startTime", "endTime"]
    opt_args = {"north":55.0 ,"south":20.0 , "east":-65.0, "west":-125.0, "coordsOrTiles":"coords", "dayNightBoth":"DNB"}

    if "-h" in args_in or "--help" in args_in or len(args_in) < 2:
        print_help(req_args, opt_args)
    elif len(args_in) % 2 != 1:
        raise IOError("All arguments must be flag-value pairs, e.g. --north 50")

    pargs = dict()
    
    for a in req_args:
        flag = "--" + a
        if flag in args_in:
            i = args_in.index(flag) + 1
            pargs[a] = args_in[i]
        else:
            raise IOError("{0} must be specified".format(a))
    
    for a in opt_args.keys():
        flag = "--" + a
        if flag in args_in:
            i = args_in.index(flag) + 1
            pargs[a] = args_in[i]
        else:
            pargs[a] = opt_args[a]

    return pargs


def print_help(req_args, opt_args):
    help_str = "\n\nUsage: python automodis.py [options]\n" \
               "\nThis function will call to MODAPS web services to find MODIS files\n" \
               "that fit the desired criteria. It works very much like submitting a\n" \
               "through LAADS Web, but can be done automatically.\n"
    opts_str_1 = "\nThere are {0} required options and {1} optional options\n".format(len(req_args), len(opt_args))
    opts_str_1b = "The API at https://ladsweb.nascom.nasa.gov/data/web_services.html\n" \
                  "will describe these options; they are the parameters for searchForFiles.\n" \
                  "All options must be given as flag-value pairs, e.g:\n" \
                  "     python automodis.py --products MYD06_L2 --north 50 --south 25\n\n"

    sys.stdout.write(help_str)
    sys.stdout.write(opts_str_1)
    sys.stdout.write(opts_str_1b)
    sys.stdout.write("Required options:")
    for arg in req_args:
        arg_str = "    --{0}\n".format(arg)
        sys.stdout.write(arg_str)

    sys.stdout.write("Optional arguments (with the default values ideal for BEHR listed):\n")
    for arg in opt_args.keys():
        arg_str = "    --{0}    (def. value: {1})\n".format(arg, opt_args[arg])
        sys.stdout.write(arg_str)

    sys.stdout.write("\n")
    exit(0)

def write_urls(urls):
    p = os.environ['MATRUNDIR']
    filename = os.path.join(p,'modis_urls.txt')
    with open(filename,'w') as f:
        for l in urls:
            f.writelines(l+"\n")

if __name__ == "__main__":
    inargs = parse_args(sys.argv)

    #url = "http://modwebsrv.modaps.eosdis.nasa.gov/axis2/services/MODAPSservices"
    url = "https://modwebsrv.modaps.eosdis.nasa.gov/axis2/services/MODAPSservices"
    server = SOAPProxy(url)
    print "Retrieving file IDs"
    attempt=0
    while True:
        try:
            fileIDs = server.searchForFiles(products=inargs["products"], startTime=inargs["startTime"], endTime=inargs["endTime"],
            north=inargs["north"], south=inargs["south"], east=inargs["east"], west=inargs["west"],
            coordsOrTiles=inargs["coordsOrTiles"], dayNightBoth=inargs["dayNightBoth"])
        except:
            if attempt > 5:
                print "More than five attempts failed to retrieve file IDs. Aborting."
                raise
            else:
                print "Retrieving file IDs failed, waiting 30 sec"
                time.sleep(30)
        else:
            break
        finally:
            attempt += 1
    
    print "fileIDs has length", len(fileIDs)
    fileIDs = ",".join(fileIDs) # returned as list, need as comma separated string
    
    attempt=0
    while True:
        try:
            fileURLs = server.getFileUrls(fileIds=fileIDs)
        except:
            if attempt > 5:
                print "More than five attempts failed to retrieve file URLs. Aborting."
                raise
            else:
                print "Retrieving file IDs failed, waiting 30 sec"
                time.sleep(30)
        else:
            break
        finally:
            attempt += 1           
    
    #fileURLs = "\n".join(fileURLs)

    # Print to shell, will be set to variable if called appropriately
    #sys.stdout.write(fileURLs)
    write_urls(fileURLs)
