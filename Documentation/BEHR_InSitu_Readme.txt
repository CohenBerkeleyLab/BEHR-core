D-AQ BEHR PRODUCT README

Josh Laughner - 11 May 2015
jlaughner@berkeley.edu

===========
1. Overview
===========

    There are two types of files for the DISCOVER-AQ BEHR product. The first
simply has the file name format "OMI_BEHR_yyyymmdd.mat" - these are the usual
BEHR files. Each one contains two variables: "Data" and "OMI" each of which are
Matlab structures. "Data" contains the native resolution OMI pixel measurements,
while the "OMI" structure contains the 0.05 x 0.05 deg gridded data.  Both are
organized by swath, so Data(1) is the first swath for that day, Data(2) is the
second, and so on. For the "OMI" structure, the entire continental US grid is
included with each swath, however as with the HDF files available online, grid
cells outside the swath will have a fill value - in this case, 0.

    The other type of file will have the file name
"OMI_BEHR_InSitu_yyyymmdd.mat". These will only contain the "Data" structure
described above, with four additional fields:
    
        InSituAMF - the air mass factor (AMF) recalculated using the observed
P3B NO2 profile. The exact method of calculation will be described below. This
field will be a NaN for any pixel without a corresponding P3B profile.

        BEHR_R_ColumnAmountNO2Trop - the tropospheric NO2 vertical column
density calculated using the InSituAMF. Again, this will have a value of NaN for
any pixels without a coinciding P3B profile. (The "_R_" indicates "reprocessed")

        ProfileCount - an integer describing the number of profiles averaged to
obtain the in situ profile used to calculate the AMF for the pixel. A value of 0
would coincide with the NaN of the previous two fields.

        InSituFlags - a currently unused field that will eventually have quality
flags about the profiles used to calculate the AMFs for each pixel.

    Although the OMI_BEHR_InSitu_yyyymmdd.mat files are only available for days
that the P3B flew, the rest of the data for that day is complete, i.e. the
regular BEHRColumnAmountNO2Trop will be correct for all pixels, not just those
with P3B profiles in them.


=================================
2. Recommended filtering criteria
=================================

    When filtering the pixels, I look for 6 criteria:
        1) The VCD must be > 0 (values < 0 are present in the OMNO2 data, but do
not make physical sense)
        2) The vcdQualityFlags field must be an even number. This indicates that
the summary bit is not set, meaning there were no significant processing issues.
        3) Cloud fraction: BEHR contains three cloud fractions: OMI geometric,
OMI radiance, and MODIS cloud fraction.
        4) The column amount should be < 1 x 10^17. Such values are expected to
indicate that the pixel has been affected by the row anomaly.
        5) The column amount must not be a NaN.
        6) Filter for row anomaly, typically by requiring that the
XTrackQualityFlags field = 0.


================================
3. Calculation of in situ values
================================

    What follows is a description of the steps taken to match P3B NO2 profiles
to relevant OMI pixels and recalculate the AMFs for those pixels.

    1) Profiles are filtered by start time: only those between 12:00 and 15:00
local standard time are used (~1.5 hrs on either side of OMI overpass). These
are NO2 profiles measured using the TD-LIF instrument.

    2) A preliminary filter is done on all pixels in a swath to remove pixels
that clearly have no overlap, by comparing boxes with edges aligned with
latitude/longitude lines around the pixels and profiles.  This is a
computationally inexpensive test that is refined next.

    3) The lat/lon of the bottom 3 km of the profile tested using the Matlab
function "inpolygon" to determine how many of those points actually fall inside
the pixel. There must be at least 20 (similar to Hains et al. JGR 2010 p.
D05301) for the pixel and profile to be considered "coincident".

    4) For each pixel associated with this profile, the profile is extended to
the BEHR pixel surface pressure and the tropopause.  If extrapolation downward
is necessary, the median of the bottom 10 NO2 measurements is taken as the
surface concentration.  The top of the profile is filled in with the nearest
WRF-Chem NO2 profile (the same WRF-Chem profiles used in BEHR).  This WRF-Chem
profile is scaled so that the top bin of the P3B profile matches the same bin in
the WRF-Chem profile.

    5) This hybrid WRF-Chem/in situ profile is then used in place of the wholly
WRF-Chem profile as the a priori in the calculation of the AMF. The scattering
weights are determined using the same parameters as the normal implementation of
BEHR (i.e. MODIS albedo and GLOBE-derived surface pressures).

 
