function [ ver_str ] = BEHR_version( )
%VER_STR = BEHR_VERSION() Returns version string for BEHR
%   BEHR version numbering follows the format n-nX where n-n is the release
%   version (NOT the collection number) of the OMNO2 product used as the
%   basis and X is a letter indicating sequential releases of BEHR. If
%   necessary, a revision number can be added as n-nXrevn, e.g. 2-1Arev1.
%   By the convention described in the change log, revision 0 is the same
%   as not including a revision number, so 2-1A and 2-1Arev0 are the same.
%   Revisions should be used when there is not a change to the core
%   algorithm, but there is some correction to the format of the output,
%   such as a changed to a fill value or something like that. In general,
%   revisions may affect how users parse the data but NOT the value of that
%   data - a user with 2-1Arev0 and 2-1Arev1 should get the same results as
%   long as they parse it properly.
%
%   Unfortunately, the release version does not seem to be well indicated
%   anywhere in the product, so we need to take this based on the last
%   publication describing updates to the product. As of 12 May 2016, that
%   is Bucsela et al. 2013, AMT, p. 2607 and the version number is 2.1. We
%   replace the . with a - so that we don't create any filename issues.
%   Whenever a substantial milestone in BEHR is reached, several changes
%   should be made:
%
%       1) This file should be updated (obviously)
%       2) The Git commit in the master branch where this version is
%       finalized should be tagged as vn.nX (e.g. v2.1A)
%       3) The data record will need to be reprocessed and published. I
%       strongly recommend that you use the cluster for this, as it can do
%       several years of BEHR in a day currently (yay parallelization).
%
%   Reference: http://disc.sci.gsfc.nasa.gov/Aura/additional/documentation/OMI_processing_info_issue1_3.pdf

ver_str = 'v2-1C';


end

