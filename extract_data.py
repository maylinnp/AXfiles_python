## Extracts data from AX titrations
# Data is saved in one structure per file
# Created 2020/09/08 MLP based on old extract-data codes
# Updated 2020/11/22 MLP: I have edited the Labview code, so that it
# applies volume correction per every 5 mL.  I also edited it so that on
# each line, the volume listed is what was added, waited 15 s, then pH on
# the corresponding line recorded (i.e., volume and pH value correspond to
# eachother). 
import csv
from itertools import zip_longest

def NaOH_calibration_data(filename: str) -> dict:
    FWD_data = list()
    # # Open filename and extract data
    with open(filename, 'r') as datafile:
        csvreader = csv.reader(datafile)
        sample_info = next(csvreader).split(",")
        print(sample_info)
        # check if FWD
        for row in csvreader:
            if "BWD" in row:
                print("There is back titration data")
                BWD_data = list()
                break
            else:
                FWD_data.append([row])
        if 'BWD_data' in locals():
            for row in csvreader:
                BWD_data.append([row])


    Lines = []
 
    w0 = sample_info[1]
    I = sample_info[2]
    emf0 = sample_info[3]
    t0 = sample_info[4]
    CNaOH = sample_info[5]
    CHCl = sample_info[6]
    if sample_info[7]:
        NaOH_ID = sample_info[7]
    else:
        NaOH_ID = None
 
    headers = ["time", "emf", "t_sample", "weight", "pH_est", "volume", "t_HCl", "t_NaOH", "t_air"]

    # process and pivot
    if FWD_data:
        FWD_data_processed = dict(zip(headers, zip_longest(*FWD_data, fillvalue=None)))

    if BWD_data:
        BWD_data_processed = dict(zip(headers, zip_longest(*BWD_data, fillvalue=None)))

    # # See if there is BWD titration data present
    Index = find(contains(Lines,'BWD')) 
    if isempty(Index)
        Index = length(Lines) 
        A = 0 
    else
        A = 1 
    end
    # # Extract FWD titration data
    clear i 
        R         = strsplit(Lines{2,1},',') 
        wHCl = str2double(R{6})/1000 

    clear R 
    # # BWD data, if it exists
    if A == 0
        emf = nan 
        t   = nan 
        w   = nan 
    else
        for i = Index+1:length(Lines)
        R         = strsplit(Lines{i,1},',') 
        emf(i-Index) = str2double(R{2}) 
        t(i-Index)   = str2double(R{3}) 
        V2uncorr(i-Index)   = str2double(R{6}) 
        Bur2T(i-Index)= str2double(R{8}) 
        end
    end
        V2 = 1.006434*V2uncorr - 0.000267*V2uncorr.**2  #dosimat 12 correction function
        rho_NaOH = 1.03121 - 1.172e-4*Bur2T - 4e-6*Bur2T.**2  
        w = V2.*rho_NaOH # NaOH #J density function 
        w = w/1000 
    end

return [w0,I,emf0,t0,CNaOH,CHCl,wHCl,emf,t,w,NaOH_ID]