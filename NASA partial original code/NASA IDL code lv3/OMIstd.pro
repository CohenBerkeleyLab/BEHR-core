
;;   New improved NO2 uncertainty calculation.  2006-11-14

PRO OmiStd, terrainReflectivity,  cloudFraction,   cloudRadianceFraction,  AMFUnpolluted,                 $  ;;Inputs
            AMFPollutedCloudy,    AMFPollutedClear,   AMFPolluted,      AMFPollutedToGround,              $
            columnAmountNO2Unpolluted, columnAmountNO2Polluted, columnAmountNO2, SlantColumnAmountNO2Std, $

            AMFUnpollutedStd, AMFPollutedStd, AMFPollutedToGroundStd, columnAmountNO2UnpollutedStd,       $  ;;Outputs
            columnAmountNO2PollutedStd,   columnAmountNO2Std,     columnAmountNO2TropStd,                 $
            columnAmountNO2BelowCloudStd, columnAmountNO2InitStd, columnAmountNO2Init, columnAmountNO2BelowCloud
;written by Eric Bucsela
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;; Computes uncertainties in OMI NO2 AMFs and vertical column amounts
;; according to the new formulas given in the document OMI_NO2_uncertainties.doc.
;; (Note, the output columnAmountNO2Init is just a de-striped version of this field).
;; These new versions are corrections to the formulas that were implemented in
;; the Public Release of OMI data in October 2006.
;;
;; Note: As of Fall 2006, uncertainties in vertical column amounts are still
;; dominated by the large values of SlantColumnAmountNO2Std.  Using the values
;; of SlantColumnAmountNO2Std given in the Level-2 file, which are on the order
;; of 2.0e15 to 2.5e15, may give unreasonably large vertical column uncertainties.
;; Based on observed statistics of SlantColumnAmountNO2, the values of
;; SlantColumnAmountNO2Std are about 1.3e15.  A similar analysis of
;; the NRT product gives slant column uncertainties on the order of 0.7e15.
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;; A priori quantities and assumptions ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  x                      = cloudFraction*0.0
  tropFractionUnpolluted = x + 0.053        ; this is approximate fixed value from the a priori unpolluted profile
  cloudAlbedo            = x + 0.8          ; cloud albedo is 80%
  AMFUnpollutedClear     = AMFUnpolluted
  AMFUnpollutedCloudy    = AMFUnpolluted




;;##############  Compute AMF uncertainties #############################################################################

;------------------------------------------------------------------------------;
; Local versions of variables that need to be scaled from values in file.

  terrainAlbedo =  (terrainReflectivity   * 1.e-3) > 0.
  cloudFrac     =  (cloudFraction         * 1.e-3) > 0.
  cloudRadFrac  =  (cloudRadianceFraction * 1.e-3) > 0.

;------------------------------------------------------------------------------;
; A priori uncertainties based on research, educated guess, or the other kind.

  terrainAlbedoStd             = x + 0.015                                     ;Koelemeijer et al.(2003)/Boersma et al.[~2007]
  cloudFracStd                 = x + 0.02                                      ;Boersma et al. [~2007]

  fracBC                       = x                                             ;fraction of polluted vertical column below cloud

  q = WHERE(cloudFrac GT 0 AND AMFPolluted GT 0) & IF MIN(q) GE 0 THEN $
  fracBC(q) = (((1. - (AMFPollutedToGround(q) / AMFPolluted(q))) / cloudFrac(q)) > 0.) < 1.

  AMFPollutedCloudyStd         = AMFPollutedCloudy * fracBC * 0.8                ;model study empirical fit

  AMFPollutedClearStd          = AMFPollutedClear / 300.   $
                               / (terrainAlbedo - 2*terrainAlbedoStd + 0.12)^2  ;model study empirical fit

  AMFUnpollutedClearStd        = AMFUnpollutedClear      * 0.02                  ;model study conservative estimate

  AMFUnpollutedCloudyStd       = AMFUnpollutedCloudy     * 0.02                  ;model study conservative estimate

  tropFractionUnpollutedStd    = tropFractionUnpolluted  * 0.20                  ;guess

  Icloud                       = cloudAlbedo   + 0.30                            ;empirical based on OMI data
  Iclear                       = terrainAlbedo + 0.25                            ;empirical based on OMI data

;------------------------------------------------------------------------------;
; Useful ratios:  ratioBC = 1 / (1 - f*fracBC) ;    ratio0 = f / w ;   ratio1  = (1-f)/(1-w)

  ratio0   = x + 1.
  ratio1   = x + 1.
  ratioBC  = x + 1.

  q = WHERE(AMFPollutedToGround GT 0)  &  IF MIN(q) GE 0 THEN $
  ratioBC  =  ((AMFPolluted / AMFPollutedToGround) > 1. ) < 1.e6    ;;ratioBC = 1/(1 - cloudFrac*fracBC) = columnPollutedGround/ColumnPolluted

  q = WHERE(cloudRadFrac LE 0. OR cloudFrac LE 0.)  & IF MIN(q) GE 0 THEN BEGIN
       ratio0(q) =  Iclear(q) / Icloud(q)
       ratio1(q) =  1.
  ENDIF

  q = WHERE(cloudRadFrac GE 1. OR cloudFrac GE 1.) & IF MIN(q) GE 0 THEN BEGIN
       ratio0(q) =  1.
       ratio1(q) =  Icloud(q) / Iclear(q)
  ENDIF

  q = WHERE(      cloudRadFrac GT 0. AND cloudRadFrac LT 1.           $
      AND cloudFrac    GT 0. AND cloudFrac    LT 1. ) & IF MIN(q) GE 0 THEN BEGIN
       ratio0(q) =    cloudFrac(q)    /    cloudRadFrac(q)
       ratio1(q) = (1.-cloudFrac(q))  /  (1.-cloudRadFrac(q))
  ENDIF


;------------------------------------------------------------------------------;
; Compute standard deviation of cloudRadianceFraction

  cloudRadFracStd = SQRT(  ( cloudRadFrac * (1.-cloudFrac) * ratio1   $
                             * terrainAlbedoStd / Icloud   )^2        $
                         +  cloudFracStd^2 ) /  (ratio1 * ratio0)


;------------------------------------------------------------------------------;
; Compute derived standard deviations of the AMF values

  AMFUnpollutedStd                  =     AMFUnpollutedCloudyStd


  AMFPollutedToGroundStd            =                                                 $
    SQRT (                                                                            $
           ( cloudRadFracStd        * (AMFPollutedCloudy-AMFPollutedClear)        )^2 $
         + ( AMFPollutedCloudyStd   * cloudRadFrac                                )^2 $
         + ( AMFPollutedClearStd    * (1. - cloudRadFrac)                         )^2 $
          )


  AMFPollutedStd                    =                                                  $
    SQRT (                                                                             $
           ( cloudRadFracStd      * (  AMFPollutedCloudy  - AMFPollutedClear           $
                                     + AMFPolluted * fracBC * ratio0 * ratio1 )   )^2  $
         + ( AMFPollutedCloudyStd * (  cloudRadFrac   -   cloudFrac * (1.-fracBC)      $
                                              * AMFPolluted / AMFPollutedCloudy ) )^2  $
         + ( AMFPollutedClearStd  * (  1. - cloudRadFrac  )                       )^2  $
         ) * ratioBC


;------------------------------------------------------------------------------;










;;##############  Compute SCD uncertainties #############################################################################


   unp = WHERE(columnAmountNO2Polluted LE 0)     ;; Identify indices of unpolluted cases

   columnAmountNO2UnpollutedStd = ABS(columnAmountNo2Unpolluted) * 0.05    ;; approx

   AMFInit = AMFUnpolluted       &       AMFInitStd = AMFUnpollutedStd

   ;; Extract columnAmountNO2Init and its uncertainty from existing de-striped fields (except SlantColumnAmountNO2Std)
   columnAmountNo2Init     = columnAmountNO2Polluted * AMFPolluted / AMFInit  +  columnAmountNo2Unpolluted
   IF MIN(unp) GE 0 THEN $
   columnAmountNo2Init(unp)= columnAmountNO2(unp)
   columnAmountNo2InitStd  = SQRT((SlantColumnAmountNO2Std / AMFInit)^2  +  (columnAmountNO2Init*AMFInitStd/AMFInit)^2)



   ;Uncertainty in polluted column amount ;;;;;;;;;;;;;;;;;;

    variance =  (columnAmountNo2UnpollutedStd^2 + columnAmountNo2InitStd^2)*   AMFInit^2 $
             +  (columnAmountNo2Unpolluted^2    - columnAmountNo2Init^2   )*AMFInitStd^2 $
             +  (columnAmountNo2Polluted        * AMFPollutedStd ) ^2

    columnAmountNo2PollutedStd = 1.-2.^15  +  x

    q = WHERE(variance GT 0 AND AMFPolluted GT 0)
    IF MIN(q) GE 0 THEN  columnAmountNo2PollutedStd(q) = SQRT(variance(q))/AMFPolluted(q)

    IF MIN(unp) GE 0 THEN $
    columnAmountNo2PollutedStd(unp) = 0.






    ;Uncertainty in tropospheric column amount ;;;;;;;;;;;;;;;;;

    variance = columnAmountNo2PollutedStd^2                                   $
             + (columnAmountNo2Unpolluted *tropFractionUnpollutedStd)^2       $
             + (tropFractionUnpolluted^2 - 2*tropFractionUnpolluted*AMFInit/AMFPolluted)*columnAmountNo2UnpollutedStd^2

    columnAmountNo2TropStd = 1.-2.^15  +  x

    q = WHERE(variance GT 0)
    IF MIN(q) GE 0 THEN columnAmountNo2TropStd(q) = SQRT(variance(q))

    IF MIN(unp) GE 0 THEN $
    columnAmountNo2TropStd(unp) =  SQRT (  (columnAmountNo2InitStd(unp) * tropFractionUnpolluted(unp))^2     $
                                         + (columnAmountNo2Init(unp)    * tropFractionUnpollutedStd(unp))^2  )






    ;Uncertainty in total column amount ;;;;;;;;;;;;;;;;;

    variance =  columnAmountNo2PollutedStd^2                                 $
             + (1 - 2*AMFInit/AMFPolluted) * columnAmountNo2UnpollutedStd^2

    columnAmountNo2Std = 1.-2.^15  +  x

    q = WHERE(variance GT 0)
    IF MIN(q) GE 0 THEN columnAmountNo2Std(q)  = SQRT(variance(q))

    IF MIN(unp) GE 0 THEN $
    columnAmountNo2Std(unp)  =  columnAmountNo2InitStd(unp)







    ;Ghost column amount and uncertainty ;;;;;;;;;;;;;;;;;;;;;

    columnAmountNo2BelowCloud    =  x
    columnAmountNo2BelowCloudStd =  x


    validGhost = (     columnAmountNo2Init                 GT columnAmountNo2Unpolluted  $
                  AND  AMFPollutedToGround                 GT            0               $
                  AND  AMFPolluted                         GT            0               $
                  AND (AMFPolluted - AMFPollutedToGround)  GT          0.001             )

    q = WHERE(validGhost)

    IF MIN(q) GE 0 THEN $
    columnAmountNo2BelowCloud(q) = (columnAmountNo2Init(q) - columnAmountNo2Unpolluted(q)) * AMFInit(q)  $
                                 * (1./AMFPollutedToGround(q) - 1./AMFPolluted(q))


    variance =    (   (columnAmountNo2UnpollutedStd^2 + columnAmountNo2InitStd^2)    * AMFInit^2        $
                    + (columnAmountNo2Unpolluted^2    - columnAmountNo2Init^2   )    * AMFInitStd^2 )   $
                * (    1./AMFPolluted                  - 1./AMFpollutedToGround   )^2                   $
             +    (    columnAmountNo2Unpolluted       - columnAmountNo2Init      )^2 * AMFInit^2       $
                * (    AMFPollutedStd/AMFPolluted^2   - AMFPollutedToGroundStd/AMFPollutedToGround^2)^2

    columnAmountNo2BelowCloudStd =  1.-2.^15  +  x

    q = WHERE(validGhost  AND  variance GE 0)
    IF MIN(q) GE 0 THEN columnAmountNo2BelowCloudStd(q) =  SQRT(variance(q))


RETURN
END
