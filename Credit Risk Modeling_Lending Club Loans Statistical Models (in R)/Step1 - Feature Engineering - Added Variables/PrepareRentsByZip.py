### PrepareRentsByZip.py
###
### Author: Bill Grieser, Team Cool French Name
###
### DESCRIPTION
###
### This script creates a csv file that contains information for each first-three-digits zip code area.
### The information includes the likely state and county for the zip code, a representative latitude and longitude,
### and 50th percentile rents for in the county for various housing sizes.
###
### To create this, the script consults the list of all US zip codes from the geonames web site. Zip codes
### in this list are grouped together based on their first three digits. Then the State abd County associated
### with the most number of zip codes in the group is assigned to the 3 digit group and a representative
### latitude and longitude are calculated. The County and State for the groups are then looked up
### in HUD 50th percentile rent data to determine the median rent in the county associated with
### the zip code.

import os
import numpy as np
import pandas as pd

### INPUT FILE CONSTANTS
### ADJUST TO FILE LOCATIONS ACCORDING TO WHERE THEY ARE ON YOUR SYSTEM
###
ZIP_CODE_FILE = "./US.txt"
MEDIAN_RENT_FILE = "./FY2012_FMRS_50_County.xlsx"

### OUTPUT FILE CONSTANTS
ZIP3_RENTS_FILE_NAME = "Zip3Rents.csv"

### OUTPUT LOCATIONS  ############################################################################################
###
### To change the working directory, set OUTPUT_FOLDER to a non-blank value that is the desired output folder
OUTPUT_FOLDER = ""
# OUTPUT_FOLDER = "D:\\Users\\billg_000\\GradSchool\\OneDrive - gwmail.gwu.edu\\DATS 6101 Intro to Data Science\\Scraping"

# Change the Current directory if a non-default output folder is selected
if OUTPUT_FOLDER != "":
    os.chdir(OUTPUT_FOLDER)

#
# Main procedure
#
def main():
    """
    Perform the main actions. This is called once the entire page is loaded.
    :return: Write files that
    """
    rents_by_zip3_df = 0

    # Create a dataframe with the rent data for Zip3 clusters
    rents_by_zip3_df = create_zip_and_rent_dataframe()

    # Write the DataFrame to a csv file, with no index column, overwriting whatever is there
    rents_by_zip3_df.to_csv(ZIP3_RENTS_FILE_NAME, index=True, mode="w")

    print ("Program ended normally.")

#
# create_zip_and_rent_dataframe
#
def create_zip_and_rent_dataframe():
    """
    Read the zip code and rent-by-county files and create a Pandas DataFrame
    with the median rent information per zip3.
    :return: A DataFrame indexed by zip3 and state
    """
    print("== Create a data frame for zip3 values and their likely location ===")

    # Define the column names for the zip code file from geonames
    col_names = ['CountryCode','ZipCode','City', 'StateName','StateCode','CountyName', 'CountyID', 'Admin3Name',
                'Admin3Code', 'Latitude', 'Longitude','LocationAcc']

    # Read the Zip code file -- tab delimited text file
    zips = pd.read_table(ZIP_CODE_FILE, sep='\t', header=None,names =col_names ,dtype="str")

    # Convert the lat and long columns to numeric
    zips["Latitude"] = pd.to_numeric(zips["Latitude"], errors="ignore")
    zips["Longitude"] = pd.to_numeric(zips["Longitude"], errors="ignore")
    zips["CountyID"] = pd.to_numeric(zips["CountyID"], errors="ignore")

    # Make a dictionary to hold the cluster information.
    # The key is the 3 digit zip code + 2 character state code, and the values are lists of zip code rows
    cluster_details = {}

    # Visit each zip. Look for zip codes that start with 0
    for zipcode in zips.iterrows():
        if len(zipcode[1].ZipCode) != 5:

            print("Illegal zip code skipped:", zipcode[1].ZipCode)
        else:
            # Get the first three characters of the zip code with the state to uniquely identify a cluster
            zip_cluster_id = zipcode[1].ZipCode[0:3] + zipcode[1].StateCode

            # See if this cluster has been encountered. If not, add an entry in
            # the cluster details
            if not zip_cluster_id in cluster_details:
                # Add a new entry for the cluster; the value is an empty list
                cluster_details[zip_cluster_id] = []

            # Add this row to the zip codes in this cluster
            cluster_details[zip_cluster_id].append(zipcode[1])

    # Make a list to hold the rows for each zip3 group
    zip3_rows = []

    # Now that the clusters have been built, visit each cluster and derive the location from it
    for zip_cluster_id in sorted(cluster_details):

        # Determine the likely location and other parameters for the zip3 group
        zip3_row = summarize_cluster(zip_cluster_id, cluster_details[zip_cluster_id])

        # Add this row to all the rows for the dataframe
        zip3_rows.append(zip3_row)

    # Make a dataframe based on the summary rows for zip codes
    zip3_df = pd.DataFrame(zip3_rows)

    print("== Reading rent data ==")

    # Visit each summary row, adding the matching rent information
    rent_df = pd.read_excel(MEDIAN_RENT_FILE)
    rent_df["county"] = pd.to_numeric(rent_df["county"])

    # Define the columns in the rent data that we care about
    # wanted_rent_cols = ["state_alpha", "county", 'Rent50_0', 'Rent50_1', 'Rent50_2', 'Rent50_3', 'Rent50_4', 'pop2000', 'countyname']
    wanted_rent_cols = ["state_alpha", "county", 'Rent50_2', 'pop2000']

    # The rent data has some counties subdivided so we need to aggregate it back to the grain we need,
    # which is one row per county
    rent_groupby = rent_df[wanted_rent_cols].groupby(["state_alpha", "county"], as_index=False)

    # Aggregate each numeric column into a new data frame
    rents_by_county = rent_groupby.median()

    # Add the county name from the rent file, just as a check
    # rents_by_county["countyname"] = rent_groupby.first()["countyname"]

    # Population in the aggregated dataframe is a sum of all the subdivisions.
    rents_by_county["pop2000"] = rent_groupby.sum()["pop2000"]

    # Define a quintile parameter for the population
    #rents_by_county["pop2000_quintile"] = pd.qcut(rents_by_county["pop2000"], [0, 0.25, 0.5, 0.75, 1],  labels=False)
    #rents_by_county["Rent50_2__quintile"] = pd.qcut(rents_by_county["Rent50_2"], [0, 0.25, 0.5, 0.75, 1],  labels=False)

    # Establish an index to use in the join to the zip3 data
    rents_by_county.set_index(["state_alpha", "county"], inplace=True)
    rents_by_county.query('state_alpha == "VA" & county == 107')

    print("== Merging Zip code and Rent-by-county data ==")

    # Join the rents into the zip file, effectively looking up the rents for each county
    # in the zips dataframe
    merged = zip3_df.join(rents_by_county, how="left", on=['StateCode', 'CountyID'])

    # Set the index on the merged dataframe. Each row pertains to one szip3/State combo
    merged.set_index(['Zip3','StateCode' ], inplace=True)

    # Return the dataframe to the caller
    return merged


### summarize_cluster
###
### Given a zip code cluster and all the zips in that cluster, return the summary information
### used to characterize that cluster. County, lat /lon
def summarize_cluster(zip_cluster, zips_in_cluster):
    """

    :param zip_cluster: The 5 character zip_cluster ID -- 3 chars for the first three digits of a zip code,
    and a two-character state code
    :param zips_in_cluster: A Python list of zip code records for all the indivual zip codes that
    start with the zip code digits and are also in the state matching the state code

    :return: A dictionary with calculated values for the zip3 group.
    """
    #
    # Split the zip_cluster ID into zip portion and state portion
    zip3 = zip_cluster[0:3]
    state_code = zip_cluster[3:]

    # Define clusters where the results are known to print out for debug
    debug_clusters = [] #['967', '606', '274', '201']

    # Debug
    if zip3 in debug_clusters:
        print(zip3, state_code, len(zips_in_cluster))

    location_counts = {}
    location_ids = {}

    'Build a list of lats and lons'
    lats = []
    lons = []

    # Visit each of the zip codes in this cluster and tally the locations
    # and acculate the list of lats and lons.
    for zip in zips_in_cluster:
        # Extract the county name from the data row
        county_name = zip.CountyName

        # If this is the first time encountering this county, we need to add
        # a starting value for this location, otherwise increment the count
        if county_name not in location_counts:
            # First time so new count = 1
            location_counts[county_name] = 1

            # Store the ID for this county to lookup later
            location_ids[county_name] = zip.CountyID

        else:
            # Increment the counts for this county
            location_counts[county_name] += 1

        # Keep a list of all the lats and lons for the whole cluster
        lats.append(zip.Latitude)
        lons.append(zip.Longitude)

    # Find the average lat and lon
    mean_lat = sum(lats) / len(lats)
    mean_lon = sum(lons) / len(lons)

    # Find the county with the max number of times referenced
    top_count = -1
    top_location = ""

    # Find the location used most often for zip codes in this cluster
    for location in location_counts:

        # See if this location is referenced more times that the current top
        if location_counts[location] > top_count:

            # If so, update the current top notion
            top_count = location_counts[location]
            top_location = location

    # Prepare a Python dictionary to return
    return_values = {}

    # Add the fields to identify this record
    return_values["StateCode"] = state_code
    return_values["Zip3"] = zip3 + "xx"

    # Add the calculated summary values to the return
    return_values['CountyName'] = top_location
    return_values['StateCode'] = state_code
    return_values['Latitude'] = mean_lat
    return_values['Longitude'] = mean_lon
    return_values['CountyID'] = location_ids[top_location]

    # Location confidence is simply the percentage of zip codes in the zip3 area that match the
    # location we are returning
    return_values['LocationConfidence'] = float(top_count) / float(len(zips_in_cluster))

    # Debug
    if zip3 in debug_clusters:
        print(location_counts)
        print(return_values)

    return return_values


# Call the main function once the file is loaded
if __name__ == "__main__":
    # Run the main program
    main()
