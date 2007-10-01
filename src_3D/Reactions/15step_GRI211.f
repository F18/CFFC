      SUBROUTINE CKWYP15STEP211(P, T, Y, WDOT)
C      SUBROUTINE CKWYP(P, T, Y, ICKWRK, RCKWRK, WDOT)
C
C  INPUT
C     P      - Pressure.
C                   cgs units - dynes/cm**2
C                   Data type - real scalar
C     T      - Temperature.
C                   cgs units - K
C                   Data type - real scalar
C     Y      - Mass fractions of the species.
C                   Data type - real array
C                   Dimension Y(*) at least 19, the total number of
C                   species.
C     Y(1)=   H2 
C     Y(2)=   H 
C     Y(3)=   O2 
C     Y(4)=   OH 
C     Y(5)=   H2O 
C     Y(6)=   HO2 
C     Y(7)=   H2O2
C     Y(8)=   CH3 
C     Y(9)=   CH4 
C     Y(10)=  CO 
C     Y(11)=  CO2 
C     Y(12)=  CH2O
C     Y(13)=  C2H2 
C     Y(14)=  C2H4
C     Y(15)=  C2H6 
C     Y(16)=  NH3 
C     Y(17)=  NO 
C     Y(18)=  HCN 
C     Y(19)=  N2
C
C     ICKWRK - Dummy Array of integer workspace.
C                   Data type - integer array 
C              (Not used - simply to be compatible with Chemkin)
C     RCKWRK - Dummy Array of real work space.
C                   Data type - real array
C              (Not used - simply to be compatible with Chemkin)
C
C  OUTPUT
C     WDOT   - Chemical molar production rates of the species.
C                   cgs units - moles/(cm**3*sec)
C                   Data type - real array
C                   Dimension WDOT(*) at least KK, the total number of
C                   species.
C
C  END PROLOGUE
C
C*****double precision
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END double precision
C*****single precision
C        IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END single precision
C
      PARAMETER (IGLO=15,IREAC=277,KK=19,KSS=29)                                
      PARAMETER (NITER=50)
      DIMENSION Y(KK), WDOT(KK)
      DIMENSION XM(IREAC), RF(IREAC), RB(IREAC), W(IREAC)
      DIMENSION RKF(IGLO), XCON(KK), WT(KK)
      DIMENSION B(KSS), ABV(KSS), DEN(KSS)
      LOGICAL LITER, LBOUND(KSS)
      SAVE RF,RB,TOLD
      DATA SMALL/1.D-50/,TOLD/0.D0/
      DATA RU/8.314510D7/
      DATA WT/  2.016,  1.008, 31.999, 17.007, 18.015, 33.007, 34.015,          
     &         15.035, 16.043, 28.011, 44.010, 30.026, 26.038, 28.054,          
     &         30.070, 17.031, 30.006, 27.026, 28.013/                          
C
C
      SUMYOW = 0.0
      DO K = 1, KK
        SUMYOW = SUMYOW + Y(K)/WT(K)
      ENDDO
      SUMYOW = SUMYOW*T*RU
      DO K = 1, KK
        XCON(K) = P*Y(K) / ( SUMYOW*WT(K) )
      ENDDO
C
C   SET LOCAL VALUES FOR THE CONCENTRATIONS
C
      XH2 = MAX( XCON(1), SMALL )                                               
      XH = MAX( XCON(2), SMALL )                                                
      XO2 = MAX( XCON(3), SMALL )                                               
      XOH = MAX( XCON(4), SMALL )                                               
      XH2O = MAX( XCON(5), SMALL )                                              
      XHO2 = MAX( XCON(6), SMALL )                                              
      XH2O2 = MAX( XCON(7), SMALL )                                             
      XCH3 = MAX( XCON(8), SMALL )                                              
      XCH4 = MAX( XCON(9), SMALL )                                              
      XCO = MAX( XCON(10), SMALL )                                              
      XCO2 = MAX( XCON(11), SMALL )                                             
      XCH2O = MAX( XCON(12), SMALL )                                            
      XC2H2 = MAX( XCON(13), SMALL )                                            
      XC2H4 = MAX( XCON(14), SMALL )                                            
      XC2H6 = MAX( XCON(15), SMALL )                                            
      XNH3 = MAX( XCON(16), SMALL )                                             
      XNO = MAX( XCON(17), SMALL )                                              
      XHCN = MAX( XCON(18), SMALL )                                             
      XN2 = MAX( XCON(19), SMALL )                                              

      BIG = 0.0
      DO 20 K = 1, KK
        IF( XCON(K) .GT. 0.0 ) THEN
          BIG = MAX( BIG, XCON(K) )
        ENDIF
        WDOT(K) = 0.0
 20   CONTINUE
C
C   EFFECTIVE THIRD BODY FOR ALL REACTIONS
C
      CALL THIRDB( XM, XH2, XH, XO2, XOH, XH2O, XHO2, XH2O2, XCH3, XCH4,        
     &             XCO, XCO2, XCH2O, XC2H2, XC2H4, XC2H6, XNH3, XNO,            
     &             XHCN, XN2 )                                                  
C
C   SET THE ELEMENTARY RATES
C
      CALL ELEMRATE( RF, RB, T, XM, TOLD )
C
C   EXPRESSIONS FOR STEADY-STATE SPECIES
C
      DO 40 K = 1, KSS
        B(K) = 0.0
        LBOUND(K) = .TRUE.
 40   CONTINUE
      CALL UPVALUE( 1, B, BIG, LBOUND, CONMAX, XNNH, XC2H, XCN, XHOCN,          
     &              XHNO, XHCCOH, XC2H5, XCH2S, XC2H3, XCH2OH, XCH2,            
     &              XHNCO, XNH2, XCH, XC, XN, XH2CN, XHCNN, XHCNO, XHCO,        
     &              XCH3OH, XCH3O, XO, XN2O, XHCCO, XNCO, XNH, XNO2,            
     &              XCH2CO )                                                    
      ADJ = 1.D0/BIG
      DO 30 N = 1, NITER

          LITER = .TRUE.
          CONMAX = 0.0
        CALL STEADYE( ABV, DEN, RF, RB, XM, ADJ, CONMAX, SMALL, LBOUND,
     &                LITER, XH2, XH, XO, XO2, XOH, XH2O, XHO2, XH2O2,          
     &                XC, XCH, XCH2, XCH2S, XCH3, XCH4, XCO, XCO2, XHCO,        
     &                XCH2O, XCH2OH, XCH3O, XCH3OH, XC2H, XC2H2, XC2H3,         
     &                XC2H4, XC2H5, XC2H6, XHCCO, XCH2CO, XHCCOH, XN,           
     &                XNH, XNH2, XNH3, XNNH, XNO, XNO2, XN2O, XHNO, XCN,        
     &                XHCN, XH2CN, XHCNN, XHCNO, XHOCN, XHNCO, XNCO,            
     &                XN2 )                                                     

        IF( LITER .AND. CONMAX .LT. 1.0D-5 ) GO TO 35


 30   CONTINUE

 35   CONTINUE
C
C   NET PRODUCTION RATES FOR SKELETAL MECHANISM
C
      CALL NETRATE( W, RF, RB, XM, XH2, XH, XO, XO2, XOH, XH2O, XHO2,           
     &              XH2O2, XC, XCH, XCH2, XCH2S, XCH3, XCH4, XCO, XCO2,         
     &              XHCO, XCH2O, XCH2OH, XCH3O, XCH3OH, XC2H, XC2H2,            
     &              XC2H3, XC2H4, XC2H5, XC2H6, XHCCO, XCH2CO, XHCCOH,          
     &              XN, XNH, XNH2, XNH3, XNNH, XNO, XNO2, XN2O, XHNO,           
     &              XCN, XHCN, XH2CN, XHCNN, XHCNO, XHOCN, XHNCO, XNCO,         
     &              XN2 )                                                       
C
C   GLOBAL RATES FOR REDUCED MECHANISM
C
      RKF(1) = +W(2) -W(9) +W(13) -W(15) +W(16) +W(19) +W(22) +W(24)            
     &         +W(25) -W(30) -W(32) +W(33) +W(34) +W(35) +W(36) +W(37)          
     &         +W(38) +W(39) +W(40) +W(41) +W(45) +W(48) + 2.00*W(49)           
     &         +W(51) + 2.00*W(53) +W(54) +W(56) -W(57) +W(58) +W(59)           
     &         +W(60) -W(63) -W(66) +W(68) +W(71) +W(72) + 3.00*W(73)           
     &         +W(75) - 2.00*W(76) +W(78) +W(80) +W(82) -W(83) -W(85)           
     &         -W(86) -W(87) -W(88) -W(89) -W(90) +W(91) - 2.00*W(95)           
     &         -W(96) -W(97) -W(98) - 2.00*W(100) -W(102) -W(103)               
     &         -W(105) +W(109) -W(111) -W(112) -W(113) + 2.00*W(116)            
     &         +W(118) -W(120) - 2.00*W(123) -W(125) +W(126)                    
     &         - 2.00*W(128) - 2.00*W(129) -W(132) +W(134) +W(135)              
     &         +W(136) -W(137) +W(138) +W(139) -W(141) +W(145) -W(146)          
     &         - 2.00*W(147) -W(149) -W(150) +W(151) +W(153) +W(154)            
     &         +W(156) +W(158) -W(159) +W(161) +W(166) +W(167) -W(170)          
     &         + 2.00*W(171) - 2.00*W(173) +W(174) -W(178) +W(181)              
     &         +W(187) +W(189) -W(190) +W(193) +W(196) -W(199) -W(201)          
     &         -W(202) -W(203) -W(206) -W(208) +W(211) +W(212) +W(214)          
     &         -W(215) - 2.00*W(216) -W(218) -W(219) -W(222) -W(229)            
     &         -W(230) +W(231) - 2.00*W(232) - 2.00*W(233) -W(234)              
     &         -W(235) +W(237) +W(238) +W(239) +W(240) +W(242) +W(244)          
     &         + 2.00*W(248) +W(249) -W(250) +W(251) -W(254) -W(255)            
     &         -W(257) - 2.00*W(258) -W(259) +W(262) +W(263) +W(264)            
     &         -W(268) +W(269) +W(272) - 2.00*W(273) +W(275) +W(277)            
      RKF(2) = +W(3) +W(9) +W(15) -W(25) +W(32) -W(38) -W(39) -W(40)            
     &         -W(41) -W(44) -W(46) -W(48) -W(49) -W(52) -W(53) -W(54)          
     &         -W(59) -W(64) -W(67) -W(68) -W(72) - 2.00*W(73) -W(74)           
     &         +W(76) -W(77) -W(78) -W(79) -W(80) +W(83) -W(91) +W(95)          
     &         +W(100) -W(109) -W(116) +W(120) +W(123) +W(125) -W(126)          
     &         +W(128) +W(129) +W(132) -W(136) +W(137) -W(138) +W(141)          
     &         +W(142) +W(143) +W(144) +W(146) + 2.00*W(147) +W(149)            
     &         +W(150) +W(159) +W(170) -W(171) -W(172) + 2.00*W(173)            
     &         -W(189) -W(195) -W(200) -W(207) -W(212) +W(215) +W(216)          
     &         +W(218) +W(219) +W(229) +W(230) +W(232) +W(233) +W(234)          
     &         +W(235) -W(237) -W(238) -W(239) -W(240) -W(242) -W(244)          
     &         -W(248) +W(250) +W(252) +W(254) +W(255) +W(257) +W(258)          
     &         +W(259) -W(264) -W(269) +W(273) -W(274) -W(275)                  
      RKF(3) = +W(4) -W(32) -W(33) -W(34) -W(35) -W(36) +W(43) +W(44)           
     &         +W(45) +W(47) -W(84) +W(86) +W(114) +W(115) +W(116)              
     &         +W(117) +W(118) +W(119) -W(166) -W(167) -W(168) -W(173)          
     &         -W(182) +W(184) -W(204) -W(214)                                  
      RKF(4) = +W(5) +W(46) +W(47) -W(84) +W(87) +W(88) -W(114) -W(115)         
     &         -W(120) +W(155)                                                  
      RKF(5) = +W(6) - 2.00*W(9) + 2.00*W(13) + 2.00*W(14) - 2.00*W(15)         
     &         +W(20) +W(21) +W(24) +W(25) +W(28) +W(29) - 2.00*W(32)           
     &         +W(48) + 2.00*W(49) -W(50) + 2.00*W(53) + 2.00*W(54)             
     &         - 2.00*W(57) -W(70) +W(72) + 3.00*W(73) - 3.00*W(76)             
     &         +W(78) +W(79) + 2.00*W(80) -W(90) + 2.00*W(91) +W(92)            
     &         - 2.00*W(95) + 2.00*W(99) - 2.00*W(100) +W(105)                  
     &         + 2.00*W(109) +W(110) +W(113) + 2.00*W(116)                      
     &         - 2.00*W(120) - 2.00*W(123) -W(124) -W(125) +W(126)              
     &         +W(127) - 2.00*W(128) - 2.00*W(129) -W(131) -W(132)              
     &         + 2.00*W(135) + 2.00*W(136) -W(137) + 2.00*W(138)                
     &         - 2.00*W(141) - 2.00*W(146) - 3.00*W(147) - 2.00*W(149)          
     &         - 2.00*W(150) + 2.00*W(158) - 2.00*W(159) + 2.00*W(164)          
     &         + 2.00*W(165) + 2.00*W(166) +W(171) +W(172)                      
     &         - 3.00*W(173) +W(174) - 2.00*W(235) +W(238) +W(239)              
     &         + 2.00*W(240) +W(244) +W(245) -W(246) + 2.00*W(247)              
     &         + 2.00*W(248) + 2.00*W(249) - 2.00*W(254) - 2.00*W(257)          
     &         - 2.00*W(258) - 2.00*W(259) +W(272) - 2.00*W(273)                
      RKF(6) = +W(7) -W(23) +W(24) +W(25) -W(30) +W(48) +W(49) -W(70)           
     &         +W(72) +W(73) -W(76) -W(89) +W(91) +W(92) -W(95) +W(110)         
     &         +W(116) -W(121) -W(123) -W(125) +W(127) -W(128) -W(129)          
     &         +W(134) +W(135) + 2.00*W(136) +W(138) +W(139) -W(141)            
     &         -W(146) -W(147) -W(149) -W(150) +W(171) +W(172) -W(173)          
     &         -W(235) -W(237) +W(240) -W(242) -W(243) +W(247) +W(248)          
     &         +W(249) -W(254) -W(259) -W(273)                                  
      RKF(7) = +W(8) +W(9) +W(15) -W(25) +W(32) -W(49) +W(50) -W(53)            
     &         +W(57) - 2.00*W(73) + 2.00*W(76) -W(78) -W(80) -W(82)            
     &         -W(91) +W(95) +W(100) -W(109) -W(116) +W(120) +W(123)            
     &         -W(126) +W(128) +W(129) +W(132) -W(135) +W(137) -W(138)          
     &         +W(141) +W(142) +W(143) +W(146) + 2.00*W(147) +W(149)            
     &         +W(150) +W(159) -W(171) + 2.00*W(173) +W(215) +W(216)            
     &         +W(218) +W(229) +W(230) +W(232) +W(233) +W(234) +W(235)          
     &         -W(237) -W(238) -W(239) -W(240) -W(242) -W(244) -W(248)          
     &         +W(250) +W(252) +W(254) +W(255) +W(257) +W(258) +W(259)          
     &         -W(269) +W(273)                                                  
      RKF(8) = +W(10) -W(15) -W(32) +W(53) -W(57) -W(60) -W(61) -W(65)          
     &         -W(66) +W(73) -W(76) +W(82) +W(91) +W(93) +W(94) -W(100)         
     &         +W(116) +W(118) -W(120) +W(126) -W(132) +W(145) +W(151)          
     &         +W(153) +W(154) +W(156) +W(157) -W(159) +W(171) -W(173)          
      RKF(9) = +W(11) +W(15) -W(16) +W(18) +W(27) +W(32) -W(49) -W(51)          
     &         +W(52) -W(53) +W(55) +W(57) -W(58) -W(59) -W(60) +W(63)          
     &         +W(66) +W(67) -W(71) -W(73) +W(74) -W(75) +W(76) +W(77)          
     &         -W(80) -W(82) -W(91) -W(93) +W(95) +W(96) +W(97) +W(100)         
     &         -W(101) +W(103) -W(109) +W(111) +W(112) -W(116) -W(117)          
     &         +W(120) +W(123) -W(126) +W(128) +W(129) +W(132) -W(135)          
     &         -W(144) -W(145) -W(151) -W(155) -W(156) -W(158) -W(161)          
     &         -W(167) -W(171) +W(172) +W(173) +W(202) +W(203) +W(204)          
     &         +W(205) +W(206) +W(207) +W(208) +W(253) +W(254) +W(273)          
     &         +W(274)                                                          
      RKF(10) = +W(12) +W(14) +W(30) +W(31) +W(98) +W(119) -W(131)              
     &          -W(151) +W(224) +W(227) +W(260) +W(266)                         
      RKF(11) = +W(26) -W(73) +W(76) -W(156) -W(157) +W(173)                    
      RKF(12) = +W(42) +W(43) +W(47) +W(61) +W(66) +W(83) +W(85) +W(86)         
     &          +W(87) +W(88) +W(92) +W(95) +W(96) +W(97) +W(99)                
     &          +W(100) +W(101) +W(102) +W(103) +W(104) +W(108) +W(110)         
     &          +W(111) +W(112) +W(113) -W(126) +W(143) -W(145) +W(191)         
     &          -W(195) +W(201) +W(208) +W(213) -W(217) +W(253) +W(265)         
     &          -W(275) -W(277)                                                 
      RKF(13) = +W(176) -W(180) -W(188) -W(189) -W(191) -W(193) +W(210)         
     &          -W(211) -W(212) -W(213) -W(214) +W(216) +W(218) -W(220)         
     &          -W(222) -W(224) -W(225) +W(229) +W(230) +W(232) +W(233)         
     &          +W(234) +W(235) -W(237) -W(238) -W(239) +W(245) +W(247)         
     &          +W(249) +W(250) +W(252) +W(254) +W(255) +W(257) +W(258)         
     &          +W(259) -W(269) +W(272) +W(273)                                 
      RKF(14) = +W(177) +W(178) -W(179) -W(181) -W(182) -W(183) -W(189)         
     &          -W(191) -W(196) +W(206) +W(216) +W(218) -W(225) -W(227)         
     &          +W(229) +W(230) +W(232) +W(233) +W(234) +W(235) -W(237)         
     &          -W(238) -W(242) -W(243) -W(244) -W(246) -W(248) -W(251)         
     &          -W(253) -W(269) +W(273)                                         
      RKF(15) = +W(198) +W(199) +W(200) +W(201) -W(234) -W(249) -W(252)         
     &          -W(263) -W(266) +W(268) +W(269) -W(272)                         
C
C   REDUCED MECHANISM
C
C  H +  0.50 O2 = OH                                                            
C  H2 +  0.50 O2 = H + OH                                                       
C  HO2 =  0.50 O2 + OH                                                          
C   0.50 O2 + H2O2 = OH + HO2                                                   
C   0.50 O2 +  0.50 C2H2 = H + CO                                               
C  CH3 + CO + C2H4 =  0.50 O2 + CH4 +  1.50 C2H2                                
C   0.50 O2 +  2CH3 = H2 + CH4 + CO                                             
C   0.50 O2 + CH3 = H + CH2O                                                    
C   0.50 O2 + CH4 = OH + CH3                                                    
C   0.50 O2 + CO = CO2                                                          
C   0.50 O2 + C2H6 = CH4 + CH2O                                                 
C  H + OH = H2O                                                                 
C  H + CH4 + NO + HCN =  0.50 O2 +  2CH3 + N2                                   
C  H +  0.50 O2 + CH4 + HCN =  2CH3 + NO                                        
C   0.50 O2 + CH4 + NH3 + HCN = H2O +  2CH3 + N2                                
C
C   SPECIES PRODUCTION RATES
C
C H2                                                                            
      WDOT(1) = - RKF(2) + RKF(7)                                               
C H                                                                             
      WDOT(2) = - RKF(1) + RKF(2) + RKF(5) + RKF(8) - RKF(12) - RKF(13)         
     &          - RKF(14)                                                       
C O2                                                                            
      WDOT(3) = - 0.50*RKF(1) - 0.50*RKF(2) + 0.50*RKF(3) - 0.50*RKF(4)         
     &          - 0.50*RKF(5) + 0.50*RKF(6) - 0.50*RKF(7) - 0.50*RKF(8)         
     &          - 0.50*RKF(9) - 0.50*RKF(10) - 0.50*RKF(11)                     
     &          + 0.50*RKF(13) - 0.50*RKF(14) - 0.50*RKF(15)                    
C OH                                                                            
      WDOT(4) = + RKF(1) + RKF(2) + RKF(3) + RKF(4) + RKF(9) - RKF(12)          
C H2O                                                                           
      WDOT(5) = + RKF(12) + RKF(15)                                             
C HO2                                                                           
      WDOT(6) = - RKF(3) + RKF(4)                                               
C H2O2                                                                          
      WDOT(7) = - RKF(4)                                                        
C CH3                                                                           
      WDOT(8) = - RKF(6) - 2.00*RKF(7) - RKF(8) + RKF(9) + 2.00*RKF(13)         
     &          + 2.00*RKF(14) + 2.00*RKF(15)                                   
C CH4                                                                           
      WDOT(9) = + RKF(6) + RKF(7) - RKF(9) + RKF(11) - RKF(13)                  
     &          - RKF(14) - RKF(15)                                             
C CO                                                                            
      WDOT(10) = + RKF(5) - RKF(6) + RKF(7) - RKF(10)                           
C CO2                                                                           
      WDOT(11) = + RKF(10)                                                      
C CH2O                                                                          
      WDOT(12) = + RKF(8) + RKF(11)                                             
C C2H2                                                                          
      WDOT(13) = - 0.50*RKF(5) + 1.50*RKF(6)                                    
C C2H4                                                                          
      WDOT(14) = - RKF(6)                                                       
C C2H6                                                                          
      WDOT(15) = - RKF(11)                                                      
C NH3                                                                           
      WDOT(16) = - RKF(15)                                                      
C NO                                                                            
      WDOT(17) = - RKF(13) + RKF(14)                                            
C HCN                                                                           
      WDOT(18) = - RKF(13) - RKF(14) - RKF(15)                                  
C N2                                                                            
      WDOT(19) = + RKF(13) + RKF(15)                                            

      DO K = 1, KK
        IF( XCON(K) .LE. 0.D0 .AND. WDOT(K) .LT. 0.D0 ) WDOT(K) = 0.D0
      END DO

      RETURN
      END
      SUBROUTINE THIRDB( XM, XH2, XH, XO2, XOH, XH2O, XHO2, XH2O2, XCH3,        
     &                   XCH4, XCO, XCO2, XCH2O, XC2H2, XC2H4, XC2H6,           
     &                   XNH3, XNO, XHCN, XN2 )                                 

C*****double precision
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END double precision
C*****single precision
C       IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END single precision

      DIMENSION XM(*)
C
C   EFFECTIVE THIRD BODY FOR ALL REACTIONS
C
      XM(1) = + 2.400*XH2+XH+XO2+XOH+15.400*XH2O+XHO2+XH2O2+XCH3                
     &        + 2.000*XCH4+ 1.750*XCO+ 3.600*XCO2+XCH2O+XC2H2+XC2H4             
     &        + 3.000*XC2H6+XNH3+XNO+XHCN+XN2                                   
      XM(2) = + 2.000*XH2+XH+XO2+XOH+ 6.000*XH2O+XHO2+XH2O2+XCH3                
     &        + 2.000*XCH4+ 1.500*XCO+ 2.000*XCO2+XCH2O+XC2H2+XC2H4             
     &        + 3.000*XC2H6+XNH3+XNO+XHCN+XN2                                   
      XM(12) = + 2.000*XH2+XH+ 6.000*XO2+XOH+ 6.000*XH2O+XHO2+XH2O2+XCH3        
     &         + 2.000*XCH4+ 1.500*XCO+ 3.500*XCO2+XCH2O+XC2H2+XC2H4            
     &         + 3.000*XC2H6+XNH3+XNO+XHCN+XN2                                  
      XM(33) = +XH2+XH+XOH+XHO2+XH2O2+XCH3+XCH4+ 0.750*XCO+ 1.500*XCO2          
     &         +XCH2O+XC2H2+XC2H4+ 1.500*XC2H6+XNH3+XNO+XHCN                    
      XM(38) = +XH+XO2+XOH+XHO2+XH2O2+XCH3+ 2.000*XCH4+XCO+XCH2O+XC2H2          
     &         +XC2H4+ 3.000*XC2H6+XNH3+XNO+XHCN+XN2                            
      XM(42) = + 0.730*XH2+XH+XO2+XOH+ 3.650*XH2O+XHO2+XH2O2+XCH3               
     &         + 2.000*XCH4+XCO+XCO2+XCH2O+XC2H2+XC2H4+ 3.000*XC2H6+XNH3        
     &         +XNO+XHCN+XN2                                                    
      XM(49) = + 2.000*XH2+XH+XO2+XOH+ 6.000*XH2O+XHO2+XH2O2+XCH3               
     &         + 2.000*XCH4+ 1.500*XCO+ 2.000*XCO2+XCH2O+XC2H2+XC2H4            
     &         + 3.000*XC2H6+XNH3+XNO+XHCN+XN2                                  
      XM(51) = + 2.000*XH2+XH+XO2+XOH+ 6.000*XH2O+XHO2+XH2O2+XCH3               
     &         + 2.000*XCH4+ 1.500*XCO+ 2.000*XCO2+XCH2O+XC2H2+XC2H4            
     &         + 3.000*XC2H6+XNH3+XNO+XHCN+XN2                                  
      XM(53) = + 2.000*XH2+XH+XO2+XOH+ 6.000*XH2O+XHO2+XH2O2+XCH3               
     &         + 2.000*XCH4+ 1.500*XCO+ 2.000*XCO2+XCH2O+XC2H2+XC2H4            
     &         + 3.000*XC2H6+XNH3+XNO+XHCN+XN2                                  
      XM(55) = + 2.000*XH2+XH+XO2+XOH+ 6.000*XH2O+XHO2+XH2O2+XCH3               
     &         + 2.000*XCH4+ 1.500*XCO+ 2.000*XCO2+XCH2O+XC2H2+XC2H4            
     &         + 3.000*XC2H6+XNH3+XNO+XHCN+XN2                                  
      XM(56) = + 2.000*XH2+XH+XO2+XOH+ 6.000*XH2O+XHO2+XH2O2+XCH3               
     &         + 2.000*XCH4+ 1.500*XCO+ 2.000*XCO2+XCH2O+XC2H2+XC2H4            
     &         + 3.000*XC2H6+XNH3+XNO+XHCN+XN2                                  
      XM(58) = + 2.000*XH2+XH+XO2+XOH+ 6.000*XH2O+XHO2+XH2O2+XCH3               
     &         + 2.000*XCH4+ 1.500*XCO+ 2.000*XCO2+XCH2O+XC2H2+XC2H4            
     &         + 3.000*XC2H6+XNH3+XNO+XHCN+XN2                                  
      XM(62) = + 2.000*XH2+XH+XO2+XOH+ 6.000*XH2O+XHO2+XH2O2+XCH3               
     &         + 2.000*XCH4+ 1.500*XCO+ 2.000*XCO2+XCH2O+XC2H2+XC2H4            
     &         + 3.000*XC2H6+XNH3+XNO+XHCN+XN2                                  
      XM(69) = + 2.000*XH2+XH+XO2+XOH+ 6.000*XH2O+XHO2+XH2O2+XCH3               
     &         + 2.000*XCH4+ 1.500*XCO+ 2.000*XCO2+XCH2O+XC2H2+XC2H4            
     &         + 3.000*XC2H6+XNH3+XNO+XHCN+XN2                                  
      XM(70) = + 2.000*XH2+XH+XO2+XOH+ 6.000*XH2O+XHO2+XH2O2+XCH3               
     &         + 2.000*XCH4+ 1.500*XCO+ 2.000*XCO2+XCH2O+XC2H2+XC2H4            
     &         + 3.000*XC2H6+XNH3+XNO+XHCN+XN2                                  
      XM(71) = + 2.000*XH2+XH+XO2+XOH+ 6.000*XH2O+XHO2+XH2O2+XCH3               
     &         + 2.000*XCH4+ 1.500*XCO+ 2.000*XCO2+XCH2O+XC2H2+XC2H4            
     &         + 3.000*XC2H6+XNH3+XNO+XHCN+XN2                                  
      XM(73) = + 2.000*XH2+XH+XO2+XOH+ 6.000*XH2O+XHO2+XH2O2+XCH3               
     &         + 2.000*XCH4+ 1.500*XCO+ 2.000*XCO2+XCH2O+XC2H2+XC2H4            
     &         + 3.000*XC2H6+XNH3+XNO+XHCN+XN2                                  
      XM(75) = + 2.000*XH2+XH+XO2+XOH+ 6.000*XH2O+XHO2+XH2O2+XCH3               
     &         + 2.000*XCH4+ 1.500*XCO+ 2.000*XCO2+XCH2O+XC2H2+XC2H4            
     &         + 3.000*XC2H6+XNH3+XNO+XHCN+XN2                                  
      XM(82) = + 2.000*XH2+XH+XO2+XOH+ 6.000*XH2O+XHO2+XH2O2+XCH3               
     &         + 2.000*XCH4+ 1.500*XCO+ 2.000*XCO2+XCH2O+XC2H2+XC2H4            
     &         + 3.000*XC2H6+XNH3+XNO+XHCN+XN2                                  
      XM(84) = + 2.000*XH2+XH+XO2+XOH+ 6.000*XH2O+XHO2+XH2O2+XCH3               
     &         + 2.000*XCH4+ 1.500*XCO+ 2.000*XCO2+XCH2O+XC2H2+XC2H4            
     &         + 3.000*XC2H6+XNH3+XNO+XHCN+XN2                                  
      XM(94) = + 2.000*XH2+XH+XO2+XOH+ 6.000*XH2O+XHO2+XH2O2+XCH3               
     &         + 2.000*XCH4+ 1.500*XCO+ 2.000*XCO2+XCH2O+XC2H2+XC2H4            
     &         + 3.000*XC2H6+XNH3+XNO+XHCN+XN2                                  
      XM(130) = + 2.000*XH2+XH+XO2+XOH+ 6.000*XH2O+XHO2+XH2O2+XCH3              
     &          + 2.000*XCH4+ 1.500*XCO+ 2.000*XCO2+XCH2O+XC2H2+XC2H4           
     &          + 3.000*XC2H6+XNH3+XNO+XHCN+XN2                                 
      XM(139) = + 2.000*XH2+XH+XO2+XOH+ 6.000*XH2O+XHO2+XH2O2+XCH3              
     &          + 2.000*XCH4+ 1.500*XCO+ 2.000*XCO2+XCH2O+XC2H2+XC2H4           
     &          + 3.000*XC2H6+XNH3+XNO+XHCN+XN2                                 
      XM(145) = + 2.000*XH2+XH+XO2+XOH+ 6.000*XH2O+XHO2+XH2O2+XCH3              
     &          + 2.000*XCH4+ 1.500*XCO+ 2.000*XCO2+XCH2O+XC2H2+XC2H4           
     &          + 3.000*XC2H6+XNH3+XNO+XHCN+XN2                                 
      XM(156) = + 2.000*XH2+XH+XO2+XOH+ 6.000*XH2O+XHO2+XH2O2+XCH3              
     &          + 2.000*XCH4+ 1.500*XCO+ 2.000*XCO2+XCH2O+XC2H2+XC2H4           
     &          + 3.000*XC2H6+XNH3+XNO+XHCN+XN2                                 
      XM(165) = + 2.000*XH2+XH+XO2+XOH+XHO2+XH2O2+XCH3+ 2.000*XCH4              
     &          + 1.500*XCO+ 2.000*XCO2+XCH2O+XC2H2+XC2H4+ 3.000*XC2H6          
     &          +XNH3+XNO+XHCN+XN2                                              
      XM(172) = + 2.000*XH2+XH+XO2+XOH+ 6.000*XH2O+XHO2+XH2O2+XCH3              
     &          + 2.000*XCH4+ 1.500*XCO+ 2.000*XCO2+XCH2O+XC2H2+XC2H4           
     &          + 3.000*XC2H6+XNH3+XNO+XHCN+XN2                                 
      XM(183) = + 2.000*XH2+XH+XO2+XOH+ 6.000*XH2O+XHO2+XH2O2+XCH3              
     &          + 2.000*XCH4+ 1.500*XCO+ 2.000*XCO2+XCH2O+XC2H2+XC2H4           
     &          + 3.000*XC2H6+XNH3+XNO+XHCN+XN2                                 
      XM(185) = + 2.000*XH2+XH+XO2+XOH+ 6.000*XH2O+XHO2+XH2O2+XCH3              
     &          + 2.000*XCH4+ 1.500*XCO+ 2.000*XCO2+XCH2O+XC2H2+XC2H4           
     &          + 3.000*XC2H6+XNH3+XNO+XHCN+XN2                                 
      XM(203) = + 2.000*XH2+XH+XO2+XOH+ 6.000*XH2O+XHO2+XH2O2+XCH3              
     &          + 2.000*XCH4+ 1.500*XCO+ 2.000*XCO2+XCH2O+XC2H2+XC2H4           
     &          + 3.000*XC2H6+XNH3+XNO+XHCN+XN2                                 
      XM(210) = + 2.000*XH2+XH+XO2+XOH+ 6.000*XH2O+XHO2+XH2O2+XCH3              
     &          + 2.000*XCH4+ 1.500*XCO+ 2.000*XCO2+XCH2O+XC2H2+XC2H4           
     &          + 3.000*XC2H6+XNH3+XNO+XHCN+XN2                                 
      XM(225) = + 2.000*XH2+XH+XO2+XOH+ 6.000*XH2O+XHO2+XH2O2+XCH3              
     &          + 2.000*XCH4+ 1.500*XCO+ 2.000*XCO2+XCH2O+XC2H2+XC2H4           
     &          + 3.000*XC2H6+XNH3+XNO+XHCN+XN2                                 
      XM(228) = + 2.000*XH2+XH+XO2+XOH+ 6.000*XH2O+XHO2+XH2O2+XCH3              
     &          + 2.000*XCH4+ 1.500*XCO+ 2.000*XCO2+XCH2O+XC2H2+XC2H4           
     &          + 3.000*XC2H6+XNH3+XNO+XHCN+XN2                                 
      XM(235) = + 2.000*XH2+XH+XO2+XOH+ 6.000*XH2O+XHO2+XH2O2+XCH3              
     &          + 2.000*XCH4+ 1.500*XCO+ 2.000*XCO2+XCH2O+XC2H2+XC2H4           
     &          + 3.000*XC2H6+XNH3+XNO+XHCN+XN2                                 
      XM(239) = + 2.000*XH2+XH+XO2+XOH+ 6.000*XH2O+XHO2+XH2O2+XCH3              
     &          + 2.000*XCH4+ 1.500*XCO+ 2.000*XCO2+XCH2O+XC2H2+XC2H4           
     &          + 3.000*XC2H6+XNH3+XNO+XHCN+XN2                                 
      XM(267) = + 2.000*XH2+XH+XO2+XOH+ 6.000*XH2O+XHO2+XH2O2+XCH3              
     &          + 2.000*XCH4+ 1.500*XCO+ 2.000*XCO2+XCH2O+XC2H2+XC2H4           
     &          + 3.000*XC2H6+XNH3+XNO+XHCN+XN2                                 

      RETURN
      END

      SUBROUTINE ELEMRATE( RF, RB, T, XM, TOLD )

C*****double precision
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END double precision
C*****single precision
C       IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END single precision

      DIMENSION RF(*), RB(*), XM(*)
      DATA SMALL/1.D-50/

C   RATE CONSTANTS FOR SKELETAL MECHANSIM

      RUC = 8.31451D0/4.184D0
      ALOGT = DLOG(T)
      RTR = 1.0D3/(RUC*T)
      TM1 = 1.0/T                                                               
      TM2 = 1.0/T/T                                                             
      TM3 = 1.0/T/T/T                                                           
      TP1 = T                                                                   
      TP2 = T*T                                                                 
      TP3 = T*T*T                                                               

      IF(ABS(T-TOLD) .GT. 1.D-3) THEN
      RF(1) = 1.2000D+17 * EXP( -1.000*ALOGT                )                   
      RB(1) = EXP(+(-3.73553D+05)*TM3+( 1.74501D+04)*TM2                        
     &        +(-5.95215D+04)*TM1+( 3.62899D+01)+(-9.65635D-04)*TP1             
     &        +( 1.49902D-07)*TP2+(-1.31833D-11)*TP3)                           
      RF(2) = 5.0000D+17 * EXP( -1.000*ALOGT                )                   
      RB(2) = EXP(+( 2.42356D+06)*TM3+(-3.77760D+04)*TM2                        
     &        +(-5.07137D+04)*TM1+( 3.39494D+01)+(-2.73498D-04)*TP1             
     &        +(-2.08792D-08)*TP2+( 4.46337D-12)*TP3)                           
      RF(3) = 5.0000D+04 * EXP(  2.670*ALOGT -  6.29068*RTR )                   
      RB(3) = EXP(+(-2.84005D+07)*TM3+( 3.35502D+05)*TM2                        
     &        +(-4.13521D+03)*TM1+( 2.85383D+01)+( 1.87759D-03)*TP1             
     &        +(-3.32578D-07)*TP2+( 3.02311D-11)*TP3)                           
      RF(4) = 2.0000D+13                                                        
      RB(4) = EXP(+(-1.27958D+06)*TM3+(-1.70702D+04)*TM2                        
     &        +(-2.65672D+04)*TM1+( 3.03272D+01)+( 3.29220D-04)*TP1             
     &        +(-5.47023D-08)*TP2+( 4.96502D-12)*TP3)                           
      RF(5) = 9.6300D+06 * EXP(  2.000*ALOGT -  4.00043*RTR )                   
      RB(5) = EXP(+(-1.75020D+07)*TM3+( 1.94410D+05)*TM2                        
     &        +(-1.05483D+04)*TM1+( 2.71317D+01)+( 2.13131D-03)*TP1             
     &        +(-3.90188D-07)*TP2+( 3.36872D-11)*TP3)                           
      RF(6) = 5.7000D+13                                                        
      RB(6) = EXP(+(-3.87771D+06)*TM3+( 5.06560D+04)*TM2                        
     &        +(-8.91233D+04)*TM1+( 3.59652D+01)+(-3.85181D-04)*TP1             
     &        +( 1.26931D-07)*TP2+(-1.29022D-11)*TP3)                           
      RF(7) = 8.0000D+13                                                        
      RB(7) = EXP(+(-7.62928D+06)*TM3+( 9.66395D+04)*TM2                        
     &        +(-4.63265D+04)*TM1+( 3.48539D+01)+(-6.64802D-04)*TP1             
     &        +( 1.64816D-07)*TP2+(-1.65701D-11)*TP3)                           
      RF(8) = 1.5000D+13                                                        
      RB(8) = EXP(+( 1.09249D+07)*TM3+(-1.20210D+05)*TM2                        
     &        +(-9.45014D+04)*TM1+( 3.22810D+01)+( 2.71096D-04)*TP1             
     &        +(-2.52275D-08)*TP2+(-3.27795D-13)*TP3)                           
      RF(9) = 1.5000D+13                                                        
      RB(9) = EXP(+( 1.28016D+06)*TM3+( 8.52154D+03)*TM2                        
     &        +(-5.05426D+04)*TM1+( 3.21800D+01)+(-5.91173D-04)*TP1             
     &        +( 1.51533D-07)*TP2+(-1.56605D-11)*TP3)                           
      RF(10) = 8.4300D+13                                                       
      RB(10) = EXP(+(-1.21447D+07)*TM3+( 1.56172D+05)*TM2                       
     &         +(-3.51673D+04)*TM1+( 3.58716D+01)+(-6.92086D-04)*TP1            
     &         +( 1.66756D-07)*TP2+(-1.61858D-11)*TP3)                          
      RF(11) = 1.0200D+09 * EXP(  1.500*ALOGT -  8.60093*RTR )                  
      RB(11) = EXP(+( 1.86553D+07)*TM3+(-2.10973D+05)*TM2                       
     &         +(-2.43174D+03)*TM1+( 2.51339D+01)+( 1.84137D-03)*TP1            
     &         +(-2.56572D-07)*TP2+( 2.02667D-11)*TP3)                          
      RF(12) = 6.0200D+14 * EXP(              -  3.00033*RTR )                  
      RB(12) = EXP(+(-7.03740D+05)*TM3+( 6.62265D+04)*TM2                       
     &         +(-6.57068D+04)*TM1+( 4.15553D+01)+(-9.34088D-04)*TP1            
     &         +( 1.42501D-07)*TP2+(-1.14781D-11)*TP3)                          
      RF(13) = 3.0000D+13                                                       
      RB(13) = EXP(+( 7.28675D+06)*TM3+(-1.04626D+05)*TM2                       
     &         +(-4.30633D+04)*TM1+( 3.04502D+01)+( 8.12574D-04)*TP1            
     &         +(-1.74620D-07)*TP2+( 1.61236D-11)*TP3)                          
      RF(14) = 3.0000D+13                                                       
      RB(14) = EXP(+( 1.39132D+07)*TM3+(-1.17252D+05)*TM2                       
     &         +(-5.58480D+04)*TM1+( 3.78861D+01)+(-5.69848D-04)*TP1            
     &         +( 1.14123D-07)*TP2+(-1.08441D-11)*TP3)                          
      RF(15) = 3.9000D+13 * EXP(              -  3.54038*RTR )                  
      RB(15) = EXP(+( 1.76497D+07)*TM3+(-2.22300D+05)*TM2                       
     &         +(-7.95215D+03)*TM1+( 2.63311D+01)+( 6.90719D-04)*TP1            
     &         +(-1.07102D-07)*TP2+( 7.61156D-12)*TP3)                          
      RF(16) = 1.0000D+13                                                       
      RB(16) = EXP(+(-1.80690D+07)*TM3+( 1.79596D+05)*TM2                       
     &         +(-3.71047D+04)*TM1+( 3.06705D+01)+( 4.61263D-04)*TP1            
     &         +(-1.02063D-07)*TP2+( 1.00951D-11)*TP3)                          
      RF(17) = 1.0000D+13                                                       
      RB(17) = EXP(+( 5.48786D+06)*TM3+(-1.06202D+05)*TM2                       
     &         +(-3.96655D+04)*TM1+( 2.71413D+01)+( 8.53658D-04)*TP1            
     &         +(-1.73931D-07)*TP2+( 1.75364D-11)*TP3)                          
      RF(18) = 3.8800D+05 * EXP(  2.500*ALOGT -  3.10034*RTR )                  
      RB(18) = EXP(+( 4.30756D+06)*TM3+(-6.40031D+04)*TM2                       
     &         +(-4.65836D+03)*TM1+( 2.46529D+01)+( 2.46874D-03)*TP1            
     &         +(-4.11511D-07)*TP2+( 3.42466D-11)*TP3)                          
      RF(19) = 1.3000D+05 * EXP(  2.500*ALOGT -  5.00054*RTR )                  
      RB(19) = EXP(+(-1.92493D+07)*TM3+( 2.21794D+05)*TM2                       
     &         +(-3.05373D+03)*TM1+( 2.70887D+01)+( 2.07635D-03)*TP1            
     &         +(-3.39643D-07)*TP2+( 2.68053D-11)*TP3)                          
      RF(20) = 5.0000D+13                                                       
      RB(20) = EXP(+(-7.92400D+06)*TM3+( 8.75558D+04)*TM2                       
     &         +(-3.97788D+04)*TM1+( 3.07605D+01)+( 7.51599D-04)*TP1            
     &         +(-1.15653D-07)*TP2+( 7.46281D-12)*TP3)                          
      RF(21) = 1.0200D+07 * EXP(  2.000*ALOGT -  1.90021*RTR )                  
      RB(21) = EXP(+(-1.39594D+07)*TM3+( 1.78089D+05)*TM2                       
     &         +(-1.20504D+04)*TM1+( 2.99686D+01)+( 1.03735D-03)*TP1            
     &         +(-1.01471D-07)*TP2+( 5.78045D-12)*TP3)                          
      RF(22) = 4.6000D+19 * EXP( -1.410*ALOGT - 28.95314*RTR )                  
      RB(22) = EXP(+( 8.50122D+06)*TM3+(-1.78101D+05)*TM2                       
     &         +( 2.14760D+03)*TM1+( 3.04179D+01)+(-2.66299D-04)*TP1            
     &         +(-1.60415D-09)*TP2+( 3.75199D-12)*TP3)                          
      RF(23) = 1.0200D+07 * EXP(  2.000*ALOGT -  1.90021*RTR )                  
      RB(23) = EXP(+(-3.62183D+07)*TM3+( 3.65800D+05)*TM2                       
     &         +(-2.61334D+04)*TM1+( 2.60321D+01)+( 2.41379D-03)*TP1            
     &         +(-4.08009D-07)*TP2+( 3.63589D-11)*TP3)                          
      RF(24) = 3.0000D+13                                                       
      RB(24) = EXP(+( 1.34855D+07)*TM3+(-1.61609D+05)*TM2                       
     &         +(-4.49385D+04)*TM1+( 3.50709D+01)+(-4.50422D-04)*TP1            
     &         +( 1.18121D-07)*TP2+(-1.22866D-11)*TP3)                          
      RF(25) = 1.9200D+07 * EXP(  1.830*ALOGT -  0.22002*RTR )                  
      RB(25) = EXP(+(-1.74594D+07)*TM3+( 1.23640D+05)*TM2                       
     &         +(-1.42993D+04)*TM1+( 2.33531D+01)+( 2.26025D-03)*TP1            
     &         +(-3.79476D-07)*TP2+( 3.25662D-11)*TP3)                          
      RF(26) = 1.3200D+14                                                       
      RB(26) = EXP(+(-8.50127D+06)*TM3+( 3.86564D+04)*TM2                       
     &         +(-3.95483D+04)*TM1+( 3.13334D+01)+( 8.34648D-04)*TP1            
     &         +(-1.55034D-07)*TP2+( 1.41614D-11)*TP3)                          
      RF(27) = 8.9800D+07 * EXP(  1.920*ALOGT -  5.69062*RTR )                  
      RB(27) = EXP(+(-2.39104D+06)*TM3+( 2.05147D+04)*TM2                       
     &         +(-4.20434D+03)*TM1+( 2.51615D+01)+( 2.16905D-03)*TP1            
     &         +(-3.70410D-07)*TP2+( 3.11994D-11)*TP3)                          
      RF(28) = 1.0000D+14                                                       
      RB(28) = EXP(+(-1.52712D+07)*TM3+( 1.00872D+05)*TM2                       
     &         +(-5.20605D+04)*TM1+( 3.04705D+01)+( 1.07588D-03)*TP1            
     &         +(-1.70100D-07)*TP2+( 1.46422D-11)*TP3)                          
      RF(29) = 1.0000D+13 * EXP(              -  8.00087*RTR )                  
      RB(29) = EXP(+( 1.24589D+07)*TM3+(-1.58599D+05)*TM2                       
     &         +(-1.48132D+03)*TM1+( 2.53621D+01)+( 7.34187D-04)*TP1            
     &         +(-1.16887D-07)*TP2+( 8.55950D-12)*TP3)                          
      RF(30) = 1.7500D+12 * EXP(              -  1.35015*RTR )                  
      RB(30) = EXP(+(-3.17357D+06)*TM3+( 1.64860D+04)*TM2                       
     &         +(-2.50023D+04)*TM1+( 2.71186D+01)+( 7.28210D-04)*TP1            
     &         +(-1.34682D-07)*TP2+( 1.21702D-11)*TP3)                          
      RF(31) = 2.5000D+12 * EXP(              - 47.80519*RTR )                  
      RB(31) = EXP(+( 9.42359D+06)*TM3+(-6.78517D+04)*TM2                       
     &         +(-2.80333D+04)*TM1+( 3.22157D+01)+(-6.90284D-04)*TP1            
     &         +( 1.17961D-07)*TP2+(-9.32107D-12)*TP3)                          
      RF(32) = 1.0000D+14 * EXP(              - 40.00434*RTR )                  
      RB(32) = EXP(+( 1.89293D+07)*TM3+(-2.05230D+05)*TM2                       
     &         +( 2.65776D+02)*TM1+( 2.75723D+01)+( 3.61499D-04)*TP1            
     &         +(-5.23996D-08)*TP2+( 2.64654D-12)*TP3)                          
      RF(33) = 2.8000D+18 * EXP( -0.860*ALOGT                )                  
      RB(33) = EXP(+( 2.33761D+06)*TM3+(-4.37788D+03)*TM2                       
     &         +(-2.42443D+04)*TM1+( 3.69366D+01)+(-5.01661D-04)*TP1            
     &         +( 1.62723D-08)*TP2+( 1.04204D-12)*TP3)                          
      RF(34) = 3.0000D+20 * EXP( -1.720*ALOGT                )                  
      RB(34) = EXP(+( 1.07259D+07)*TM3+(-1.04678D+05)*TM2                       
     &         +(-2.36435D+04)*TM1+( 3.56837D+01)+(-1.12244D-03)*TP1            
     &         +( 1.24084D-07)*TP2+(-8.44059D-12)*TP3)                          
      RF(35) = 9.3800D+18 * EXP( -0.760*ALOGT                )                  
      RB(35) = EXP(+( 1.36223D+06)*TM3+( 7.28493D+03)*TM2                       
     &         +(-2.43141D+04)*TM1+( 3.88348D+01)+(-4.29478D-04)*TP1            
     &         +( 3.73602D-09)*TP2+( 2.14467D-12)*TP3)                          
      RF(36) = 3.7500D+20 * EXP( -1.720*ALOGT                )                  
      RB(36) = EXP(+( 1.07259D+07)*TM3+(-1.04678D+05)*TM2                       
     &         +(-2.36435D+04)*TM1+( 3.59068D+01)+(-1.12244D-03)*TP1            
     &         +( 1.24084D-07)*TP2+(-8.44059D-12)*TP3)                          
      RF(37) = 8.3000D+13 * EXP(              - 14.41456*RTR )                  
      RB(37) = EXP(+( 2.79711D+06)*TM3+(-5.52261D+04)*TM2                       
     &         +( 1.55414D+03)*TM1+( 2.82823D+01)+( 6.92137D-04)*TP1            
     &         +(-1.70781D-07)*TP2+( 1.76467D-11)*TP3)                          
      RF(38) = 1.0000D+18 * EXP( -1.000*ALOGT                )                  
      RB(38) = EXP(+( 4.78151D+06)*TM3+(-6.18811D+04)*TM2                       
     &         +(-5.16092D+04)*TM1+( 3.53257D+01)+(-2.23802D-04)*TP1            
     &         +(-2.30197D-08)*TP2+( 3.67249D-12)*TP3)                          
      RF(39) = 9.0000D+16 * EXP( -0.600*ALOGT                )                  
      RB(39) = EXP(+( 8.79998D+05)*TM3+(-1.52298D+04)*TM2                       
     &         +(-5.18886D+04)*TM1+( 3.56745D+01)+( 6.49300D-05)*TP1            
     &         +(-7.31647D-08)*TP2+( 8.08301D-12)*TP3)                          
      RF(40) = 6.0000D+19 * EXP( -1.250*ALOGT                )                  
      RB(40) = EXP(+( 7.21995D+06)*TM3+(-9.10381D+04)*TM2                       
     &         +(-5.14346D+04)*TM1+( 3.76970D+01)+(-4.04260D-04)*TP1            
     &         +( 8.32098D-09)*TP2+( 9.15911D-13)*TP3)                          
      RF(41) = 5.5000D+20 * EXP( -2.000*ALOGT                )                  
      RB(41) = EXP(+( 1.45353D+07)*TM3+(-1.78509D+05)*TM2                       
     &         +(-5.09106D+04)*TM1+( 3.47436D+01)+(-9.45633D-04)*TP1            
     &         +( 1.02343D-07)*TP2+(-7.35382D-12)*TP3)                          
      RF(42) = 2.2000D+22 * EXP( -2.000*ALOGT                )                  
      RB(42) = EXP(+( 3.77110D+06)*TM3+(-4.81691D+04)*TM2                       
     &         +(-5.90404D+04)*TM1+( 4.05243D+01)+(-1.15350D-03)*TP1            
     &         +( 1.14110D-07)*TP2+(-6.69238D-12)*TP3)                          
      RF(43) = 3.9700D+12 * EXP(              -  0.67107*RTR )                  
      RB(43) = EXP(+(-6.88870D+06)*TM3+( 3.39387D+04)*TM2                       
     &         +(-2.71224D+04)*TM1+( 2.77176D+01)+( 8.63184D-04)*TP1            
     &         +(-2.15857D-07)*TP2+( 2.24823D-11)*TP3)                          
      RF(44) = 2.8000D+13 * EXP(              -  1.06812*RTR )                  
      RB(44) = EXP(+( 1.07837D+06)*TM3+(-4.11753D+04)*TM2                       
     &         +(-2.80002D+04)*TM1+( 3.13468D+01)+( 3.78915D-04)*TP1            
     &         +(-5.68427D-08)*TP2+( 4.17414D-12)*TP3)                          
      RF(45) = 1.3400D+14 * EXP(              -  0.63507*RTR )                  
      RB(45) = EXP(+( 1.51753D+06)*TM3+(-7.22963D+04)*TM2                       
     &         +(-1.80790D+04)*TM1+( 2.84618D+01)+( 1.02136D-03)*TP1            
     &         +(-2.25484D-07)*TP2+( 2.26117D-11)*TP3)                          
      RF(46) = 1.2100D+07 * EXP(  2.000*ALOGT -  5.20056*RTR )                  
      RB(46) = EXP(+(-1.51440D+07)*TM3+( 1.70305D+05)*TM2                       
     &         +(-1.20478D+04)*TM1+( 2.80431D+01)+( 2.18100D-03)*TP1            
     &         +(-3.92328D-07)*TP2+( 3.28963D-11)*TP3)                          
      RF(47) = 1.0000D+13 * EXP(              -  3.60039*RTR )                  
      RB(47) = EXP(+(-4.88311D+06)*TM3+(-4.90772D+03)*TM2                       
     &         +(-3.57346D+04)*TM1+( 2.59087D+01)+( 1.55083D-03)*TP1            
     &         +(-3.55319D-07)*TP2+( 3.41168D-11)*TP3)                          
      RF(48) = 1.1000D+14                                                       
      RB(48) = EXP(+( 5.62116D+06)*TM3+(-5.30958D+04)*TM2                       
     &         +(-1.16826D+04)*TM1+( 3.31389D+01)+( 1.42451D-04)*TP1            
     &         +( 1.39284D-08)*TP2+(-3.94506D-12)*TP3)                          
      RF(50) = 3.0000D+13                                                       
      RB(50) = EXP(+( 1.48026D+07)*TM3+(-1.70866D+05)*TM2                       
     &         +(-5.37812D+03)*TM1+( 2.86830D+01)+( 6.56277D-04)*TP1            
     &         +(-1.52159D-07)*TP2+( 1.25744D-11)*TP3)                          
      RF(52) = 6.6000D+08 * EXP(  1.620*ALOGT - 10.84118*RTR )                  
      RB(52) = EXP(+( 1.98428D+07)*TM3+(-2.21082D+05)*TM2                       
     &         +(-4.53842D+03)*TM1+( 2.62087D+01)+( 1.97768D-03)*TP1            
     &         +(-2.73756D-07)*TP2+( 2.07990D-11)*TP3)                          
      RF(54) = 7.3400D+13                                                       
      RB(54) = EXP(+( 9.64470D+06)*TM3+(-1.28731D+05)*TM2                       
     &         +(-4.39588D+04)*TM1+( 3.20280D+01)+( 8.62269D-04)*TP1            
     &         +(-1.76760D-07)*TP2+( 1.53327D-11)*TP3)                          
      RF(57) = 2.3000D+10 * EXP(  1.050*ALOGT -  3.27536*RTR )                  
      RB(57) = EXP(+( 9.76617D+06)*TM3+(-1.23946D+05)*TM2                       
     &         +(-9.44779D+03)*TM1+( 2.68150D+01)+( 1.49834D-03)*TP1            
     &         +(-2.40873D-07)*TP2+( 1.83983D-11)*TP3)                          
      RF(59) = 2.0000D+13                                                       
      RB(59) = EXP(+(-1.57111D+07)*TM3+( 1.55491D+05)*TM2                       
     &         +(-3.80002D+04)*TM1+( 3.20467D+01)+( 5.10958D-04)*TP1            
     &         +(-1.04203D-07)*TP2+( 9.30426D-12)*TP3)                          
      RF(60) = 1.2000D+13                                                       
      RB(60) = EXP(+(-5.92436D+06)*TM3+( 2.34233D+04)*TM2                       
     &         +(-1.93737D+03)*TM1+( 2.70466D+01)+( 1.15335D-03)*TP1            
     &         +(-2.68819D-07)*TP2+( 2.62810D-11)*TP3)                          
      RF(61) = 6.0000D+12                                                       
      RB(61) = EXP(+(-1.01057D+07)*TM3+( 5.50086D+04)*TM2                       
     &         +(-1.75797D+03)*TM1+( 2.61302D+01)+( 1.58498D-03)*TP1            
     &         +(-3.51071D-07)*TP2+( 3.32378D-11)*TP3)                          
      RF(63) = 3.4000D+06 * EXP(  1.600*ALOGT                )                  
      RB(63) = EXP(+( 7.95084D+06)*TM3+(-9.91922D+04)*TM2                       
     &         +(-3.67854D+03)*TM1+( 2.25372D+01)+( 1.54732D-03)*TP1            
     &         +(-2.72448D-07)*TP2+( 2.50834D-11)*TP3)                          
      RF(64) = 2.0000D+13                                                       
      RB(64) = EXP(+( 7.84581D+06)*TM3+(-1.30307D+05)*TM2                       
     &         +(-4.05610D+04)*TM1+( 2.85175D+01)+( 9.03353D-04)*TP1            
     &         +(-1.76072D-07)*TP2+( 1.67455D-11)*TP3)                          
      RF(65) = 3.2000D+13                                                       
      RB(65) = EXP(+( 1.76325D+07)*TM3+(-2.62374D+05)*TM2                       
     &         +(-4.49821D+03)*TM1+( 2.44982D+01)+( 1.54574D-03)*TP1            
     &         +(-3.40688D-07)*TP2+( 3.37223D-11)*TP3)                          
      RF(66) = 1.6000D+13                                                       
      RB(66) = EXP(+( 1.34511D+07)*TM3+(-2.30789D+05)*TM2                       
     &         +(-4.31882D+03)*TM1+( 2.35818D+01)+( 1.97738D-03)*TP1            
     &         +(-4.22939D-07)*TP2+( 4.06791D-11)*TP3)                          
      RF(67) = 1.7000D+07 * EXP(  2.100*ALOGT -  4.87053*RTR )                  
      RB(67) = EXP(+( 1.05670D+07)*TM3+(-1.34759D+05)*TM2                       
     &         +(-6.16524D+03)*TM1+( 2.63592D+01)+( 2.22970D-03)*TP1            
     &         +(-3.63506D-07)*TP2+( 2.90452D-11)*TP3)                          
      RF(68) = 4.2000D+06 * EXP(  2.100*ALOGT -  4.87053*RTR )                  
      RB(68) = EXP(+(-1.29899D+07)*TM3+( 1.51038D+05)*TM2                       
     &         +(-3.60440D+03)*TM1+( 2.84902D+01)+( 1.83731D-03)*TP1            
     &         +(-2.91638D-07)*TP2+( 2.16039D-11)*TP3)                          
      RF(72) = 3.0000D+13                                                       
      RB(72) = EXP(+( 2.27542D+07)*TM3+(-2.89145D+05)*TM2                       
     &         +(-3.35921D+04)*TM1+( 3.11358D+01)+( 7.39776D-04)*TP1            
     &         +(-1.50161D-07)*TP2+( 1.17542D-11)*TP3)                          
      RF(74) = 1.3250D+06 * EXP(  2.530*ALOGT - 12.24133*RTR )                  
      RB(74) = EXP(+(-5.19346D+06)*TM3+( 4.62776D+04)*TM2                       
     &         +(-3.33116D+03)*TM1+( 2.60762D+01)+( 2.63102D-03)*TP1            
     &         +(-4.43848D-07)*TP2+( 3.67619D-11)*TP3)                          
      RF(76) = 2.0000D+12                                                       
      RB(76) = EXP(+( 1.11164D+07)*TM3+(-1.17959D+05)*TM2                       
     &         +(-3.37042D+04)*TM1+( 2.88930D+01)+( 6.35767D-04)*TP1            
     &         +(-1.14214D-07)*TP2+( 8.59401D-12)*TP3)                          
      RF(77) = 1.1500D+08 * EXP(  1.900*ALOGT -  7.53082*RTR )                  
      RB(77) = EXP(+( 1.61981D+05)*TM3+(-5.92296D+03)*TM2                       
     &         +(-6.01191D+03)*TM1+( 2.59541D+01)+( 2.20431D-03)*TP1            
     &         +(-3.70043D-07)*TP2+( 3.01880D-11)*TP3)                          
      RF(78) = 1.0000D+14                                                       
      RB(78) = EXP(+(-3.11684D+07)*TM3+( 2.75828D+05)*TM2                       
     &         +(-9.86681D+03)*TM1+( 2.92997D+01)+( 1.30281D-03)*TP1            
     &         +(-2.93255D-07)*TP2+( 2.96688D-11)*TP3)                          
      RF(79) = 5.0000D+13 * EXP(              -  8.00087*RTR )                  
      RB(79) = EXP(+( 1.48168D+07)*TM3+(-1.82704D+05)*TM2                       
     &         +(-2.37684D+03)*TM1+( 2.76546D+01)+( 7.83882D-04)*TP1            
     &         +(-1.19028D-07)*TP2+( 7.76862D-12)*TP3)                          
      RF(80) = 1.1300D+13 * EXP(              -  3.42837*RTR )                  
      RB(80) = EXP(+(-2.29343D+07)*TM3+( 1.91879D+05)*TM2                       
     &         +(-1.82519D+04)*TM1+( 2.55460D+01)+( 1.44720D-03)*TP1            
     &         +(-3.18264D-07)*TP2+( 3.11420D-11)*TP3)                          
      RF(81) = 1.0000D+13                                                       
      RB(81) = EXP(+(-1.15433D+07)*TM3+( 1.46248D+05)*TM2                       
     &         +(-1.57368D+04)*TM1+( 2.98964D+01)+( 4.78691D-05)*TP1            
     &         +( 1.31837D-08)*TP2+(-1.74954D-12)*TP3)                          
      RF(83) = 2.1600D+08 * EXP(  1.510*ALOGT -  3.43037*RTR )                  
      RB(83) = EXP(+(-2.54924D+07)*TM3+( 3.06449D+05)*TM2                       
     &         +(-1.09109D+04)*TM1+( 3.16895D+01)+( 8.82096D-04)*TP1            
     &         +(-1.77530D-07)*TP2+( 1.73112D-11)*TP3)                          
      RF(85) = 3.5700D+04 * EXP(  2.400*ALOGT +  2.11023*RTR )                  
      RB(85) = EXP(+(-3.18153D+07)*TM3+( 3.86143D+05)*TM2                       
     &         +(-9.63997D+03)*TM1+( 2.97985D+01)+( 1.57422D-03)*TP1            
     &         +(-2.91243D-07)*TP2+( 2.63337D-11)*TP3)                          
      RF(86) = 2.9000D+13 * EXP(              +  0.50005*RTR )                  
      RB(86) = EXP(+(-9.68582D+06)*TM3+( 8.91648D+04)*TM2                       
     &         +(-3.53409D+04)*TM1+( 3.34736D+01)+( 1.71047D-04)*TP1            
     &         +(-4.50753D-08)*TP2+( 4.83558D-12)*TP3)                          
      RF(87) = 1.7500D+12 * EXP(              -  0.32003*RTR )                  
      RB(87) = EXP(+(-6.40064D+06)*TM3+( 6.73886D+04)*TM2                       
     &         +(-1.63245D+04)*TM1+( 2.82328D+01)+( 5.29473D-04)*TP1            
     &         +(-1.29836D-07)*TP2+( 1.15051D-11)*TP3)                          
      RF(88) = 5.8000D+14 * EXP(              -  9.56104*RTR )                  
      RB(88) = EXP(+(-6.40064D+06)*TM3+( 6.73886D+04)*TM2                       
     &         +(-2.09747D+04)*TM1+( 3.40362D+01)+( 5.29473D-04)*TP1            
     &         +(-1.29836D-07)*TP2+( 1.15051D-11)*TP3)                          
      RF(89) = 5.0000D+13                                                       
      RB(89) = EXP(+(-7.14092D+06)*TM3+( 7.96468D+04)*TM2                       
     &         +(-7.83362D+04)*TM1+( 3.57099D+01)+(-4.77936D-04)*TP1            
     &         +( 1.10862D-07)*TP2+(-9.74799D-12)*TP3)                          
      RF(90) = 3.0000D+13                                                       
      RB(90) = EXP(+(-1.11645D+07)*TM3+( 1.55282D+05)*TM2                       
     &         +(-4.60600D+04)*TM1+( 3.59054D+01)+(-1.19775D-03)*TP1            
     &         +( 3.01551D-07)*TP2+(-2.90258D-11)*TP3)                          
      RF(91) = 2.0000D+13                                                       
      RB(91) = EXP(+(-2.52790D+07)*TM3+( 3.18940D+05)*TM2                       
     &         +(-4.01559D+04)*TM1+( 3.84311D+01)+(-1.35552D-03)*TP1            
     &         +( 2.71918D-07)*TP2+(-2.41817D-11)*TP3)                          
      RF(92) = 1.1300D+07 * EXP(  2.000*ALOGT -  3.00033*RTR )                  
      RB(92) = EXP(+(-2.43786D+07)*TM3+( 2.80849D+05)*TM2                       
     &         +(-1.21987D+04)*TM1+( 3.07668D+01)+( 1.81844D-03)*TP1            
     &         +(-3.77833D-07)*TP2+( 3.43788D-11)*TP3)                          
      RF(93) = 3.0000D+13                                                       
      RB(93) = EXP(+(-1.63695D+07)*TM3+( 2.30822D+05)*TM2                       
     &         +(-4.43720D+04)*TM1+( 3.78366D+01)+(-1.28189D-03)*TP1            
     &         +( 2.58635D-07)*TP2+(-2.32721D-11)*TP3)                          
      RF(95) = 5.6000D+07 * EXP(  1.600*ALOGT -  5.42059*RTR )                  
      RB(95) = EXP(+(-1.08780D+07)*TM3+( 1.30072D+05)*TM2                       
     &         +(-7.88219D+03)*TM1+( 2.76448D+01)+( 1.66019D-03)*TP1            
     &         +(-2.96115D-07)*TP2+( 2.55085D-11)*TP3)                          
      RF(96) = 2.5010D+13                                                       
      RB(96) = EXP(+(-4.18137D+06)*TM3+( 3.15853D+04)*TM2                       
     &         +( 1.79394D+02)*TM1+( 3.06270D+01)+( 4.31632D-04)*TP1            
     &         +(-8.22514D-08)*TP2+( 6.95683D-12)*TP3)                          
      RF(97) = 1.0000D+08 * EXP(  1.600*ALOGT -  3.12034*RTR )                  
      RB(97) = EXP(+( 9.27365D+06)*TM3+(-9.30748D+04)*TM2                       
     &         +(-8.76901D+03)*TM1+( 2.62756D+01)+( 1.75538D-03)*TP1            
     &         +(-2.59481D-07)*TP2+( 2.12399D-11)*TP3)                          
      RF(98) = 4.7600D+07 * EXP(  1.228*ALOGT -  0.07001*RTR )                  
      RB(98) = EXP(+(-5.35116D+06)*TM3+( 1.30594D+05)*TM2                       
     &         +(-1.36778D+04)*TM1+( 3.35776D+01)+(-4.96013D-04)*TP1            
     &         +( 1.34797D-07)*TP2+(-1.34274D-11)*TP3)                          
      RF(99) = 5.0000D+13                                                       
      RB(99) = EXP(+(-1.11948D+06)*TM3+( 1.60898D+03)*TM2                       
     &         +(-5.20886D+04)*TM1+( 3.37359D+01)+( 6.54400D-04)*TP1            
     &         +(-1.64993D-07)*TP2+( 1.59942D-11)*TP3)                          
      RF(100) = 3.4300D+09 * EXP(  1.180*ALOGT +  0.44705*RTR )                 
      RB(100) = EXP(+(-2.26600D+06)*TM3+( 2.15558D+04)*TM2                      
     &          +(-1.57952D+04)*TM1+( 2.78998D+01)+( 1.38431D-03)*TP1           
     &          +(-2.45403D-07)*TP2+( 2.04932D-11)*TP3)                         
      RF(101) = 5.0000D+12                                                      
      RB(101) = EXP(+(-2.64752D+07)*TM3+( 2.85831D+05)*TM2                      
     &          +(-4.61300D+04)*TM1+( 3.27522D+01)+( 3.03090D-04)*TP1           
     &          +(-9.24361D-08)*TP2+( 9.96570D-12)*TP3)                         
      RF(102) = 5.0000D+12                                                      
      RB(102) = EXP(+(-2.91837D+06)*TM3+( 3.33709D+01)*TM2                      
     &          +(-4.86909D+04)*TM1+( 2.92230D+01)+( 6.95485D-04)*TP1           
     &          +(-1.64304D-07)*TP2+( 1.74070D-11)*TP3)                         
      RF(103) = 1.4400D+06 * EXP(  2.000*ALOGT +  0.84009*RTR )                 
      RB(103) = EXP(+( 7.78212D+05)*TM3+(-1.60821D+04)*TM2                      
     &          +(-1.13515D+04)*TM1+( 2.52932D+01)+( 1.94965D-03)*TP1           
     &          +(-3.39203D-07)*TP2+( 2.86040D-11)*TP3)                         
      RF(104) = 6.3000D+06 * EXP(  2.000*ALOGT -  1.50016*RTR )                 
      RB(104) = EXP(+(-2.27787D+07)*TM3+( 2.69715D+05)*TM2                      
     &          +(-9.96833D+03)*TM1+( 3.02983D+01)+( 1.55726D-03)*TP1           
     &          +(-2.67334D-07)*TP2+( 2.11627D-11)*TP3)                         
      RF(105) = 2.0000D+13                                                      
      RB(105) = EXP(+( 1.07997D+07)*TM3+(-4.15119D+04)*TM2                      
     &          +(-2.54294D+04)*TM1+( 3.58130D+01)+(-1.15780D-03)*TP1           
     &          +( 3.27620D-07)*TP2+(-3.55713D-11)*TP3)                         
      RF(106) = 2.1800D-04 * EXP(  4.500*ALOGT +  1.00011*RTR )                 
      RB(106) = EXP(+(-5.08027D+07)*TM3+( 6.28259D+05)*TM2                      
     &          +(-1.48822D+04)*TM1+( 2.72011D+01)+( 2.10774D-03)*TP1           
     &          +(-2.97990D-07)*TP2+( 2.47867D-11)*TP3)                         
      RF(107) = 5.0400D+05 * EXP(  2.300*ALOGT - 13.50146*RTR )                 
      RB(107) = EXP(+(-1.78011D+07)*TM3+( 2.25429D+05)*TM2                      
     &          +(-4.90595D+03)*TM1+( 3.36373D+01)+( 4.71839D-04)*TP1           
     &          +(-3.53759D-08)*TP2+( 2.27839D-12)*TP3)                         
      RF(108) = 3.3700D+07 * EXP(  2.000*ALOGT - 14.00152*RTR )                 
      RB(108) = EXP(+(-3.31654D+07)*TM3+( 3.25836D+05)*TM2                      
     &          +(-1.73592D+03)*TM1+( 2.87523D+01)+( 2.03697D-03)*TP1           
     &          +(-4.19464D-07)*TP2+( 4.12223D-11)*TP3)                         
      RF(109) = 4.8300D-04 * EXP(  4.000*ALOGT +  2.00022*RTR )                 
      RB(109) = EXP(+(-6.88602D+07)*TM3+( 7.61824D+05)*TM2                      
     &          +(-3.05563D+04)*TM1+( 2.00408D+01)+( 3.19402D-03)*TP1           
     &          +(-5.53572D-07)*TP2+( 5.04156D-11)*TP3)                         
      RF(110) = 5.0000D+12                                                      
      RB(110) = EXP(+( 1.19900D+07)*TM3+(-1.58805D+05)*TM2                      
     &          +(-4.17219D+04)*TM1+( 3.14358D+01)+( 5.31907D-04)*TP1           
     &          +(-1.38393D-07)*TP2+( 1.24157D-11)*TP3)                         
      RF(111) = 3.6000D+06 * EXP(  2.000*ALOGT -  2.50027*RTR )                 
      RB(111) = EXP(+(-1.07881D+07)*TM3+( 1.14805D+05)*TM2                      
     &          +(-6.18888D+03)*TM1+( 2.55148D+01)+( 2.04058D-03)*TP1           
     &          +(-3.65638D-07)*TP2+( 3.15794D-11)*TP3)                         
      RF(112) = 3.5400D+06 * EXP(  2.120*ALOGT -  0.87009*RTR )                 
      RB(112) = EXP(+(-1.27480D+07)*TM3+( 1.50075D+05)*TM2                      
     &          +(-1.09436D+04)*TM1+( 2.60813D+01)+( 2.15524D-03)*TP1           
     &          +(-3.85855D-07)*TP2+( 3.32753D-11)*TP3)                         
      RF(113) = 7.5000D+12 * EXP(              -  2.00022*RTR )                 
      RB(113) = EXP(+( 4.05265D+06)*TM3+(-5.23637D+04)*TM2                      
     &          +(-7.48703D+03)*TM1+( 2.78493D+01)+( 5.76014D-04)*TP1           
     &          +(-1.07260D-07)*TP2+( 8.43006D-12)*TP3)                         
      RF(114) = 1.3000D+11 * EXP(              +  1.63018*RTR )                 
      RB(114) = EXP(+(-3.28517D+06)*TM3+( 2.17763D+04)*TM2                      
     &          +(-1.86088D+04)*TM1+( 2.80240D+01)+(-3.58426D-04)*TP1           
     &          +( 8.47605D-08)*TP2+(-6.66952D-12)*TP3)                         
      RF(115) = 4.2000D+14 * EXP(              - 12.00130*RTR )                 
      RB(115) = EXP(+(-3.28517D+06)*TM3+( 2.17763D+04)*TM2                      
     &          +(-2.54683D+04)*TM1+( 3.61044D+01)+(-3.58426D-04)*TP1           
     &          +( 8.47605D-08)*TP2+(-6.66952D-12)*TP3)                         
      RF(116) = 2.0000D+13                                                      
      RB(116) = EXP(+(-2.37614D+07)*TM3+( 2.46644D+05)*TM2                      
     &          +(-5.79153D+04)*TM1+( 3.43640D+01)+(-3.34163D-04)*TP1           
     &          +( 4.64342D-08)*TP2+(-1.57001D-12)*TP3)                         
      RF(117) = 1.0000D+12                                                      
      RB(117) = EXP(+(-3.45655D+07)*TM3+( 3.68845D+05)*TM2                      
     &          +(-2.95114D+04)*TM1+( 3.32786D+01)+(-4.29403D-04)*TP1           
     &          +( 1.38259D-08)*TP2+( 1.23780D-12)*TP3)                         
      RF(118) = 2.0000D+13                                                      
      RB(118) = EXP(+(-1.61150D+07)*TM3+( 1.90078D+05)*TM2                      
     &          +(-1.32612D+04)*TM1+( 3.31582D+01)+(-5.24386D-04)*TP1           
     &          +( 1.15204D-07)*TP2+(-1.11106D-11)*TP3)                         
      RF(119) = 1.5000D+14 * EXP(              - 23.60256*RTR )                 
      RB(119) = EXP(+( 8.14400D+06)*TM3+(-8.49219D+04)*TM2                      
     &          +(-4.24214D+04)*TM1+( 3.60105D+01)+(-3.61064D-04)*TP1           
     &          +( 6.32591D-08)*TP2+(-4.35605D-12)*TP3)                         
      RF(120) = 1.0000D+12 * EXP(              -  8.00087*RTR )                 
      RB(120) = EXP(+( 1.56441D+07)*TM3+(-1.83454D+05)*TM2                      
     &          +(-3.05863D+03)*TM1+( 2.54003D+01)+( 3.07240D-06)*TP1           
     &          +( 3.23609D-08)*TP2+(-4.02297D-12)*TP3)                         
      RF(121) = 5.8000D+13 * EXP(              -  0.57606*RTR )                 
      RB(121) = EXP(+(-4.34381D+06)*TM3+( 2.44207D+04)*TM2                      
     &          +(-6.98183D+04)*TM1+( 3.20908D+01)+( 2.14201D-04)*TP1           
     &          +(-5.99191D-08)*TP2+( 7.89869D-12)*TP3)                         
      RF(122) = 5.0000D+13                                                      
      RB(122) = EXP(+( 4.31826D+06)*TM3+(-6.65517D+04)*TM2                      
     &          +(-3.88238D+04)*TM1+( 3.44601D+01)+(-6.96583D-04)*TP1           
     &          +( 8.97804D-08)*TP2+(-4.75515D-12)*TP3)                         
      RF(123) = 5.0000D+13                                                      
      RB(123) = EXP(+( 2.27042D+07)*TM3+(-2.15664D+05)*TM2                      
     &          +(-4.95676D+04)*TM1+( 3.56016D+01)+(-7.84630D-04)*TP1           
     &          +( 1.62984D-07)*TP2+(-1.60584D-11)*TP3)                         
      RF(124) = 3.3000D+13                                                      
      RB(124) = EXP(+(-8.36735D+06)*TM3+( 1.00056D+05)*TM2                      
     &          +(-3.72522D+04)*TM1+( 3.22332D+01)+(-5.05617D-04)*TP1           
     &          +( 1.30770D-07)*TP2+(-1.13791D-11)*TP3)                         
      RF(125) = 1.1070D+08 * EXP(  1.790*ALOGT -  1.67018*RTR )                 
      RB(125) = EXP(+(-2.33524D+07)*TM3+( 2.91512D+05)*TM2                      
     &          +(-9.28919D+02)*TM1+( 3.22083D+01)+( 7.09430D-04)*TP1           
     &          +(-8.55237D-08)*TP2+( 8.07232D-12)*TP3)                         
      RF(126) = 1.7130D+13 * EXP(              +  0.75508*RTR )                 
      RB(126) = EXP(+(-2.04079D+07)*TM3+( 2.71348D+05)*TM2                      
     &          +(-3.04841D+04)*TM1+( 3.75337D+01)+(-1.73030D-03)*TP1           
     &          +( 3.99026D-07)*TP2+(-3.65079D-11)*TP3)                         
      RF(127) = 4.0000D+13                                                      
      RB(127) = EXP(+( 1.28331D+07)*TM3+(-8.18875D+04)*TM2                      
     &          +(-6.53433D+04)*TM1+( 3.95008D+01)+(-1.35531D-03)*TP1           
     &          +( 2.84215D-07)*TP2+(-2.72084D-11)*TP3)                         
      RF(128) = 3.0000D+13                                                      
      RB(128) = EXP(+( 5.57118D+06)*TM3+( 2.03853D+04)*TM2                      
     &          +(-2.76581D+04)*TM1+( 3.57946D+01)+(-1.38196D-03)*TP1           
     &          +( 3.27073D-07)*TP2+(-3.17576D-11)*TP3)                         
      RF(129) = 6.0000D+13                                                      
      RB(129) = EXP(+( 2.17315D+07)*TM3+(-1.40843D+05)*TM2                      
     &          +(-3.02056D+04)*TM1+( 3.66811D+01)+(-1.37843D-03)*TP1           
     &          +( 3.83085D-07)*TP2+(-3.76867D-11)*TP3)                         
      RF(131) = 3.4000D+12 * EXP(              -  0.69007*RTR )                 
      RB(131) = EXP(+(-1.77909D+07)*TM3+( 1.67908D+05)*TM2                      
     &          +(-3.36225D+04)*TM1+( 2.62921D+01)+( 1.84667D-04)*TP1           
     &          +( 1.28082D-08)*TP2+(-2.05805D-12)*TP3)                         
      RF(132) = 9.4600D+13 * EXP(              +  0.51506*RTR )                 
      RB(132) = EXP(+( 3.12013D+07)*TM3+(-2.97396D+05)*TM2                      
     &          +(-3.71701D+04)*TM1+( 3.71755D+01)+(-1.14029D-03)*TP1           
     &          +( 2.78439D-07)*TP2+(-2.78584D-11)*TP3)                         
      RF(133) = 5.0000D+13                                                      
      RB(133) = EXP(+(-9.42584D+06)*TM3+( 1.05823D+05)*TM2                      
     &          +(-7.94262D+04)*TM1+( 3.57875D+01)+( 2.11343D-05)*TP1           
     &          +(-2.23232D-08)*TP2+( 3.37000D-12)*TP3)                         
      RF(134) = 1.3200D+13 * EXP(              -  1.50016*RTR )                 
      RB(134) = EXP(+(-4.83217D+06)*TM3+( 4.14134D+04)*TM2                      
     &          +(-3.82736D+04)*TM1+( 2.92846D+01)+( 2.73353D-05)*TP1           
     &          +(-5.96537D-09)*TP2+( 1.07654D-12)*TP3)                         
      RF(135) = 5.0000D+05 * EXP(  2.000*ALOGT -  7.23078*RTR )                 
      RB(135) = EXP(+(-3.49998D+07)*TM3+( 4.20129D+05)*TM2                      
     &          +(-9.12882D+03)*TM1+( 3.02213D+01)+( 7.30532D-04)*TP1           
     &          +(-1.43423D-07)*TP2+( 1.48476D-11)*TP3)                         
      RF(136) = 3.2000D+13                                                      
      RB(136) = EXP(+( 1.87262D+07)*TM3+(-1.64635D+05)*TM2                      
     &          +(-6.65052D+04)*TM1+( 3.79284D+01)+(-7.72662D-04)*TP1           
     &          +( 1.45339D-07)*TP2+(-1.55436D-11)*TP3)                         
      RF(137) = 4.0000D+13                                                      
      RB(137) = EXP(+(-8.01929D+06)*TM3+( 1.86429D+05)*TM2                      
     &          +(-3.34163D+04)*TM1+( 4.01904D+01)+(-1.60410D-03)*TP1           
     &          +( 3.14879D-07)*TP2+(-2.89582D-11)*TP3)                         
      RF(138) = 2.4600D+06 * EXP(  2.000*ALOGT -  8.27090*RTR )                 
      RB(138) = EXP(+( 6.44065D+05)*TM3+( 1.01091D+04)*TM2                      
     &          +(-7.60352D+03)*TM1+( 2.65506D+01)+( 1.53885D-03)*TP1           
     &          +(-2.14092D-07)*TP2+( 1.77840D-11)*TP3)                         
      RF(140) = 3.0000D+13                                                      
      RB(140) = EXP(+(-2.62869D+07)*TM3+( 3.12221D+05)*TM2                      
     &          +(-4.69961D+04)*TM1+( 3.38238D+01)+(-1.35994D-04)*TP1           
     &          +(-1.10378D-08)*TP2+( 3.28056D-12)*TP3)                         
      RF(141) = 1.5000D+13 * EXP(              -  0.60007*RTR )                 
      RB(141) = EXP(+( 8.90944D+06)*TM3+(-8.81179D+04)*TM2                      
     &          +(-4.51812D+03)*TM1+( 2.93391D+01)+( 7.36294D-05)*TP1           
     &          +(-1.32833D-08)*TP2+( 9.09603D-13)*TP3)                         
      RF(142) = 2.8000D+13                                                      
      RB(142) = EXP(+( 1.86942D+07)*TM3+(-2.30183D+05)*TM2                      
     &          +(-3.33859D+04)*TM1+( 2.83665D+01)+( 4.65205D-04)*TP1           
     &          +(-4.76266D-08)*TP2+( 2.62009D-12)*TP3)                         
      RF(143) = 1.2000D+13                                                      
      RB(143) = EXP(+( 2.95779D+06)*TM3+(-4.50956D+04)*TM2                      
     &          +(-9.38234D+04)*TM1+( 3.03821D+01)+( 7.55365D-04)*TP1           
     &          +(-1.84241D-07)*TP2+( 1.79803D-11)*TP3)                         
      RF(144) = 7.0000D+13                                                      
      RB(144) = EXP(+(-6.58281D+06)*TM3+( 9.87548D+04)*TM2                      
     &          +(-8.30921D+03)*TM1+( 3.41946D+01)+(-6.39501D-04)*TP1           
     &          +( 9.40188D-08)*TP2+(-6.29539D-12)*TP3)                         
      RF(146) = 3.0000D+13                                                      
      RB(146) = EXP(+( 8.90944D+06)*TM3+(-8.81179D+04)*TM2                      
     &          +(-4.21616D+03)*TM1+( 3.00323D+01)+( 7.36294D-05)*TP1           
     &          +(-1.32833D-08)*TP2+( 9.09603D-13)*TP3)                         
      RF(147) = 1.2000D+13 * EXP(              +  0.57006*RTR )                 
      RB(147) = EXP(+( 8.90157D+05)*TM3+( 9.83113D+04)*TM2                      
     &          +(-3.73456D+04)*TM1+( 3.79865D+01)+(-1.53047D-03)*TP1           
     &          +( 3.01595D-07)*TP2+(-2.80486D-11)*TP3)                         
      RF(148) = 1.6000D+13 * EXP(              +  0.57006*RTR )                 
      RB(148) = EXP(+( 2.90611D+07)*TM3+(-3.11265D+05)*TM2                      
     &          +(-5.97364D+03)*TM1+( 2.74546D+01)+( 1.68818D-04)*TP1           
     &          +( 2.33502D-08)*TP2+(-3.35905D-12)*TP3)                         
      RF(149) = 9.0000D+12                                                      
      RB(149) = EXP(+( 8.90944D+06)*TM3+(-8.81179D+04)*TM2                      
     &          +(-4.21616D+03)*TM1+( 2.88283D+01)+( 7.36294D-05)*TP1           
     &          +(-1.32833D-08)*TP2+( 9.09603D-13)*TP3)                         
      RF(150) = 7.0000D+12                                                      
      RB(150) = EXP(+( 8.90944D+06)*TM3+(-8.81179D+04)*TM2                      
     &          +(-4.21616D+03)*TM1+( 2.85770D+01)+( 7.36294D-05)*TP1           
     &          +(-1.32833D-08)*TP2+( 9.09603D-13)*TP3)                         
      RF(151) = 1.4000D+13                                                      
      RB(151) = EXP(+(-2.29960D+07)*TM3+( 2.43448D+05)*TM2                      
     &          +(-3.15873D+04)*TM1+( 2.96385D+01)+( 1.00530D-04)*TP1           
     &          +(-3.01081D-08)*TP2+( 3.69565D-12)*TP3)                         
      RF(152) = 4.0000D+13 * EXP(              +  0.55006*RTR )                 
      RB(152) = EXP(+( 1.21113D+07)*TM3+(-1.28762D+05)*TM2                      
     &          +(-8.92743D+03)*TM1+( 2.79339D+01)+( 1.93331D-04)*TP1           
     &          +(-3.78352D-08)*TP2+( 2.94265D-12)*TP3)                         
      RF(153) = 2.6750D+13 * EXP(              - 28.80312*RTR )                 
      RB(153) = EXP(+(-1.48354D+07)*TM3+( 2.07148D+05)*TM2                      
     &          +(-1.18821D+03)*TM1+( 3.37485D+01)+(-8.53606D-04)*TP1           
     &          +( 1.69906D-07)*TP2+(-1.60756D-11)*TP3)                         
      RF(154) = 3.6000D+10 * EXP(              -  8.94097*RTR )                 
      RB(154) = EXP(+(-9.34754D+06)*TM3+( 1.00946D+05)*TM2                      
     &          +(-3.08588D+04)*TM1+( 2.43455D+01)+( 5.15675D-08)*TP1           
     &          +(-4.02516D-09)*TP2+( 1.46085D-12)*TP3)                         
      RF(155) = 2.4500D+04 * EXP(  2.470*ALOGT -  5.18056*RTR )                 
      RB(155) = EXP(+(-5.53721D+07)*TM3+( 6.35140D+05)*TM2                      
     &          +(-1.44147D+04)*TM1+( 3.03441D+01)+( 1.71195D-03)*TP1           
     &          +(-3.80580D-07)*TP2+( 3.51423D-11)*TP3)                         
      RF(157) = 4.9900D+12 * EXP(  0.100*ALOGT - 10.60115*RTR )                 
      RB(157) = EXP(+(-4.61876D+06)*TM3+( 1.29179D+05)*TM2                      
     &          +(-1.02355D+03)*TM1+( 3.49143D+01)+(-1.45455D-03)*TP1           
     &          +( 3.09254D-07)*TP2+(-2.92446D-11)*TP3)                         
      RF(158) = 2.6480D+13                                                      
      RB(158) = EXP(+(-2.59992D+07)*TM3+( 2.81289D+05)*TM2                      
     &          +(-4.60075D+04)*TM1+( 3.62725D+01)+( 5.39505D-05)*TP1           
     &          +(-1.06092D-07)*TP2+( 1.23964D-11)*TP3)                         
      RF(159) = 3.3200D+03 * EXP(  2.810*ALOGT -  5.86064*RTR )                 
      RB(159) = EXP(+(-4.30443D+07)*TM3+( 4.91340D+05)*TM2                      
     &          +(-1.40269D+04)*TM1+( 2.84579D+01)+( 1.96044D-03)*TP1           
     &          +(-3.90843D-07)*TP2+( 3.48683D-11)*TP3)                         
      RF(160) = 3.0000D+07 * EXP(  1.500*ALOGT -  9.94108*RTR )                 
      RB(160) = EXP(+(-1.92246D+07)*TM3+( 2.05284D+05)*TM2                      
     &          +(-1.03464D+04)*TM1+( 2.80560D+01)+( 9.88287D-04)*TP1           
     &          +(-2.17620D-07)*TP2+( 1.94931D-11)*TP3)                         
      RF(161) = 1.0000D+07 * EXP(  1.500*ALOGT -  9.94108*RTR )                 
      RB(161) = EXP(+(-4.27815D+07)*TM3+( 4.91081D+05)*TM2                      
     &          +(-7.78556D+03)*TM1+( 3.04866D+01)+( 5.95892D-04)*TP1           
     &          +(-1.45752D-07)*TP2+( 1.20518D-11)*TP3)                         
      RF(162) = 2.2700D+05 * EXP(  2.000*ALOGT -  9.20100*RTR )                 
      RB(162) = EXP(+(-3.56678D+07)*TM3+( 3.94485D+05)*TM2                      
     &          +(-3.47969D+03)*TM1+( 2.59233D+01)+( 1.44013D-03)*TP1           
     &          +(-3.06737D-07)*TP2+( 2.79817D-11)*TP3)                         
      RF(163) = 6.1400D+06 * EXP(  1.740*ALOGT - 10.45113*RTR )                 
      RB(163) = EXP(+(-3.39213D+07)*TM3+( 3.85436D+05)*TM2                      
     &          +(-9.41840D+03)*TM1+( 2.71854D+01)+( 1.28050D-03)*TP1           
     &          +(-2.79316D-07)*TP2+( 2.54875D-11)*TP3)                         
      RF(164) = 2.2440D+18 * EXP( -1.000*ALOGT - 17.00184*RTR )                 
      RB(164) = EXP(+( 2.43707D+07)*TM3+(-3.00106D+05)*TM2                      
     &          +( 4.91900D+02)*TM1+( 3.46927D+01)+(-3.57591D-04)*TP1           
     &          +( 9.69846D-08)*TP2+(-1.03924D-11)*TP3)                         
      RF(165) = 1.8700D+17 * EXP( -1.000*ALOGT - 17.00184*RTR )                 
      RB(165) = EXP(+( 2.43707D+07)*TM3+(-3.00106D+05)*TM2                      
     &          +( 4.91900D+02)*TM1+( 3.22078D+01)+(-3.57591D-04)*TP1           
     &          +( 9.69846D-08)*TP2+(-1.03924D-11)*TP3)                         
      RF(166) = 7.6000D+12 * EXP(              -  0.40004*RTR )                 
      RB(166) = EXP(+( 8.56633D+06)*TM3+(-8.75559D+04)*TM2                      
     &          +(-1.66974D+04)*TM1+( 2.93767D+01)+( 4.83354D-04)*TP1           
     &          +(-1.19918D-07)*TP2+( 1.11586D-11)*TP3)                         
      RF(167) = 1.8000D+13 * EXP(              -  0.90010*RTR )                 
      RB(167) = EXP(+(-1.67894D+07)*TM3+( 1.96666D+05)*TM2                      
     &          +(-1.09904D+04)*TM1+( 3.15578D+01)+( 1.32043D-04)*TP1           
     &          +(-4.73608D-08)*TP2+( 5.13012D-12)*TP3)                         
      RF(168) = 4.2800D-13 * EXP(  7.600*ALOGT +  3.53038*RTR )                 
      RB(168) = EXP(+(-6.73612D+07)*TM3+( 7.97243D+05)*TM2                      
     &          +(-1.66308D+04)*TM1+( 2.14066D+01)+( 6.01035D-03)*TP1           
     &          +(-1.07198D-06)*TP2+( 9.63714D-11)*TP3)                         
      RF(169) = 5.0000D+13 * EXP(              -  1.50016*RTR )                 
      RB(169) = EXP(+(-1.62913D+07)*TM3+( 1.87612D+05)*TM2                      
     &          +(-7.77859D+04)*TM1+( 3.18661D+01)+( 2.45982D-04)*TP1           
     &          +( 1.51166D-08)*TP2+(-3.91631D-12)*TP3)                         
      RF(170) = 4.0700D+05 * EXP(  2.400*ALOGT -  0.20002*RTR )                 
      RB(170) = EXP(+(-2.05154D+07)*TM3+( 3.17668D+05)*TM2                      
     &          +(-1.66140D+04)*TM1+( 3.39138D+01)+( 9.31217D-04)*TP1           
     &          +(-1.20364D-07)*TP2+( 7.95494D-12)*TP3)                         
      RF(171) = 3.9800D+12 * EXP(              +  0.24003*RTR )                 
      RB(171) = EXP(+(-2.60832D+07)*TM3+( 2.35843D+05)*TM2                      
     &          +(-4.46407D+04)*TM1+( 2.91618D+01)+( 1.84252D-04)*TP1           
     &          +(-2.95475D-08)*TP2+( 4.19268D-12)*TP3)                         
      RF(173) = 8.4000D+11 * EXP(              -  3.87542*RTR )                 
      RB(173) = EXP(+( 1.00380D+07)*TM3+(-7.67841D+04)*TM2                      
     &          +(-8.19167D+03)*TM1+( 2.76420D+01)+( 2.56852D-04)*TP1           
     &          +(-5.73709D-08)*TP2+( 4.41987D-12)*TP3)                         
      RF(174) = 1.6000D+12 * EXP(              -  0.85409*RTR )                 
      RB(174) = EXP(+(-1.24741D+07)*TM3+( 4.56455D+04)*TM2                      
     &          +(-4.36825D+04)*TM1+( 2.25678D+01)+( 1.76802D-03)*TP1           
     &          +(-3.40881D-07)*TP2+( 3.22889D-11)*TP3)                         
      RF(175) = 1.0000D+13                                                      
      RB(175) = EXP(+(-2.08194D+07)*TM3+( 1.56039D+05)*TM2                      
     &          +(-4.23634D+04)*TM1+( 2.81212D+01)+( 1.48220D-03)*TP1           
     &          +(-3.19354D-07)*TP2+( 3.09144D-11)*TP3)                         
      RF(176) = 3.5000D+13 * EXP(              -  0.33004*RTR )                 
      RB(176) = EXP(+(-1.32159D+06)*TM3+( 1.39970D+03)*TM2                      
     &          +(-3.79787D+04)*TM1+( 3.23987D+01)+( 1.82461D-04)*TP1           
     &          +(-4.79577D-08)*TP2+( 5.15276D-12)*TP3)                         
      RF(177) = 2.6500D+12 * EXP(              -  6.40069*RTR )                 
      RB(177) = EXP(+(-2.41773D+06)*TM3+( 8.44426D+03)*TM2                      
     &          +(-1.90949D+04)*TM1+( 2.68618D+01)+( 1.11305D-04)*TP1           
     &          +(-2.94596D-08)*TP2+( 3.94173D-12)*TP3)                         
      RF(178) = 7.3330D+13 * EXP(              -  1.12012*RTR )                 
      RB(178) = EXP(+(-5.21484D+06)*TM3+( 6.36704D+04)*TM2                      
     &          +(-2.52454D+04)*TM1+( 3.39498D+01)+(-5.80832D-04)*TP1           
     &          +( 1.41322D-07)*TP2+(-1.37050D-11)*TP3)                         
      RF(179) = 1.4000D+12 * EXP(              - 10.81117*RTR )                 
      RB(179) = EXP(+(-1.23698D+07)*TM3+( 9.77729D+04)*TM2                      
     &          +(-4.54062D+04)*TM1+( 2.58400D+01)+( 8.25092D-04)*TP1           
     &          +(-1.52200D-07)*TP2+( 1.34169D-11)*TP3)                         
      RF(180) = 2.9000D+13 * EXP(              - 23.15251*RTR )                 
      RB(180) = EXP(+(-1.34660D+07)*TM3+( 1.04817D+05)*TM2                      
     &          +(-2.96779D+04)*TM1+( 2.59147D+01)+( 7.53937D-04)*TP1           
     &          +(-1.33702D-07)*TP2+( 1.22059D-11)*TP3)                         
      RF(181) = 4.4000D+14 * EXP(              - 18.88205*RTR )                 
      RB(181) = EXP(+(-9.57274D+06)*TM3+( 4.25468D+04)*TM2                      
     &          +(-4.06598D+04)*TM1+( 2.78227D+01)+( 1.51723D-03)*TP1           
     &          +(-3.22981D-07)*TP2+( 3.10636D-11)*TP3)                         
      RF(182) = 2.0000D+12 * EXP(              - 21.06228*RTR )                 
      RB(182) = EXP(+(-1.10903D+07)*TM3+( 1.14843D+05)*TM2                      
     &          +(-2.39976D+04)*TM1+( 2.64962D+01)+( 4.95872D-04)*TP1           
     &          +(-9.74974D-08)*TP2+( 8.45189D-12)*TP3)                         
      RF(184) = 2.1100D+12 * EXP(              +  0.48005*RTR )                 
      RB(184) = EXP(+( 6.81130D+06)*TM3+(-6.97113D+04)*TM2                      
     &          +(-3.17529D+03)*TM1+( 3.01293D+01)+(-2.19440D-04)*TP1           
     &          +( 6.05302D-08)*TP2+(-6.31833D-12)*TP3)                         
      RF(185) = 1.0600D+20 * EXP( -1.410*ALOGT                )                 
      RB(185) = EXP(+( 1.17164D+07)*TM3+(-8.30085D+04)*TM2                      
     &          +(-3.60847D+04)*TM1+( 4.22990D+01)+(-1.81025D-03)*TP1           
     &          +( 3.16533D-07)*TP2+(-2.89874D-11)*TP3)                         
      RF(186) = 3.9000D+12 * EXP(              +  0.24003*RTR )                 
      RB(186) = EXP(+(-8.09089D+06)*TM3+( 5.26411D+04)*TM2                      
     &          +(-2.30296D+04)*TM1+( 2.69409D+01)+( 5.48660D-04)*TP1           
     &          +(-1.15233D-07)*TP2+( 1.12833D-11)*TP3)                         
      RF(187) = 1.3200D+14 * EXP(              -  0.36004*RTR )                 
      RB(187) = EXP(+(-5.29377D+06)*TM3+(-2.58504D+03)*TM2                      
     &          +(-1.45237D+04)*TM1+( 2.66952D+01)+( 1.24080D-03)*TP1           
     &          +(-2.86014D-07)*TP2+( 2.89300D-11)*TP3)                         
      RF(188) = 5.0000D+13                                                      
      RB(188) = EXP(+(-2.71384D+06)*TM3+( 4.19597D+04)*TM2                      
     &          +(-3.59236D+04)*TM1+( 3.40763D+01)+(-4.70927D-04)*TP1           
     &          +( 1.18676D-07)*TP2+(-1.13940D-11)*TP3)                         
      RF(189) = 3.2000D+13 * EXP(              -  0.33004*RTR )                 
      RB(189) = EXP(+( 4.85895D+06)*TM3+(-4.58158D+04)*TM2                      
     &          +(-1.23035D+04)*TM1+( 3.22893D+01)+( 1.59601D-04)*TP1           
     &          +(-2.47862D-08)*TP2+( 1.52004D-12)*TP3)                         
      RF(190) = 2.0000D+13                                                      
      RB(190) = EXP(+(-8.93667D+06)*TM3+( 1.33256D+05)*TM2                      
     &          +(-9.31154D+03)*TM1+( 3.53906D+01)+(-7.47463D-04)*TP1           
     &          +( 9.53621D-08)*TP2+(-6.88167D-12)*TP3)                         
      RF(191) = 2.0000D+09 * EXP(  1.200*ALOGT                )                 
      RB(191) = EXP(+(-1.76098D+07)*TM3+( 2.24478D+05)*TM2                      
     &          +(-2.11055D+04)*TM1+( 3.29711D+01)+( 8.17929D-04)*TP1           
     &          +(-1.63454D-07)*TP2+( 1.54131D-11)*TP3)                         
      RF(192) = 4.6100D+05 * EXP(  2.000*ALOGT -  6.50071*RTR )                 
      RB(192) = EXP(+(-2.56471D+07)*TM3+( 3.11287D+05)*TM2                      
     &          +(-5.17214D+03)*TM1+( 2.78214D+01)+( 1.38834D-03)*TP1           
     &          +(-3.26144D-07)*TP2+( 3.28176D-11)*TP3)                         
      RF(193) = 1.2800D+06 * EXP(  1.500*ALOGT -  0.10001*RTR )                 
      RB(193) = EXP(+(-1.45474D+07)*TM3+( 1.61676D+05)*TM2                      
     &          +(-2.82140D+04)*TM1+( 2.31660D+01)+( 1.30396D-03)*TP1           
     &          +(-2.40149D-07)*TP2+( 2.27921D-11)*TP3)                         
      RF(194) = 1.5000D+13                                                      
      RB(194) = EXP(+(-4.03543D+06)*TM3+( 4.33594D+04)*TM2                      
     &          +(-7.37363D+04)*TM1+( 3.40846D+01)+(-2.88467D-04)*TP1           
     &          +( 7.07184D-08)*TP2+(-6.24127D-12)*TP3)                         
      RF(195) = 2.0000D+13 * EXP(              - 13.85150*RTR )                 
      RB(195) = EXP(+( 1.82752D+06)*TM3+( 2.91633D+03)*TM2                      
     &          +(-8.15203D+03)*TM1+( 3.32988D+01)+(-5.39594D-04)*TP1           
     &          +( 8.35947D-08)*TP2+(-7.54311D-12)*TP3)                         
      RF(196) = 2.1600D+13 * EXP( -0.230*ALOGT                )                 
      RB(196) = EXP(+( 3.42278D+06)*TM3+(-4.71355D+04)*TM2                      
     &          +(-4.88938D+04)*TM1+( 3.08403D+01)+( 1.26345D-04)*TP1           
     &          +(-4.17700D-08)*TP2+( 4.92763D-12)*TP3)                         
      RF(197) = 4.1600D+14 * EXP( -0.450*ALOGT                )                 
      RB(197) = EXP(+( 1.51414D+07)*TM3+(-1.15340D+05)*TM2                      
     &          +(-1.75821D+04)*TM1+( 3.81772D+01)+(-1.54969D-03)*TP1           
     &          +( 3.08791D-07)*TP2+(-2.85617D-11)*TP3)                         
      RF(198) = 7.0000D+12                                                      
      RB(198) = EXP(+( 1.00627D+07)*TM3+(-1.21559D+05)*TM2                      
     &          +(-4.89776D+03)*TM1+( 2.76738D+01)+( 3.87491D-04)*TP1           
     &          +(-6.02043D-08)*TP2+( 4.86765D-12)*TP3)                         
      RF(199) = 4.6000D+13                                                      
      RB(199) = EXP(+( 1.12603D+06)*TM3+( 1.16972D+04)*TM2                      
     &          +(-1.42093D+04)*TM1+( 3.43204D+01)+(-3.59971D-04)*TP1           
     &          +( 3.51579D-08)*TP2+(-2.01402D-12)*TP3)                         
      RF(200) = 4.0000D+13 * EXP(              -  3.65040*RTR )                 
      RB(200) = EXP(+( 1.24206D+07)*TM3+(-1.45664D+05)*TM2                      
     &          +(-7.63022D+03)*TM1+( 3.00998D+01)+( 4.37186D-04)*TP1           
     &          +(-6.23447D-08)*TP2+( 4.07677D-12)*TP3)                         
      RF(201) = 9.0000D+07 * EXP(  1.500*ALOGT +  0.46005*RTR )                 
      RB(201) = EXP(+(-1.29742D+07)*TM3+( 1.59618D+05)*TM2                      
     &          +(-1.47394D+04)*TM1+( 2.95250D+01)+( 1.31206D-03)*TP1           
     &          +(-2.38621D-07)*TP2+( 2.12777D-11)*TP3)                         
      RF(202) = 3.3000D+08                                                      
      RB(202) = EXP(+( 7.78579D+06)*TM3+(-1.26913D+05)*TM2                      
     &          +(-3.47247D+03)*TM1+( 1.98525D+01)+( 3.60247D-04)*TP1           
     &          +(-3.43906D-08)*TP2+( 1.85266D-12)*TP3)                         
      RF(203) = 1.3000D+14 * EXP( -0.110*ALOGT -  4.98054*RTR )                 
      RB(203) = EXP(+( 8.85871D+06)*TM3+(-1.39742D+05)*TM2                      
     &          +(-5.90192D+03)*TM1+( 3.19784D+01)+( 2.80846D-04)*TP1           
     &          +(-2.06007D-08)*TP2+( 6.39763D-13)*TP3)                         
      RF(204) = 5.0000D+12                                                      
      RB(204) = EXP(+( 1.73516D+06)*TM3+(-3.09903D+04)*TM2                      
     &          +(-2.83175D+04)*TM1+( 2.98660D+01)+( 4.79361D-04)*TP1           
     &          +(-1.25930D-07)*TP2+( 1.23773D-11)*TP3)                         
      RF(205) = 2.5000D+13                                                      
      RB(205) = EXP(+( 4.55577D+05)*TM3+(-4.80605D+04)*TM2                      
     &          +(-5.48847D+04)*TM1+( 3.11759D+01)+( 8.08581D-04)*TP1           
     &          +(-1.80632D-07)*TP2+( 1.73423D-11)*TP3)                         
      RF(206) = 7.0000D+13                                                      
      RB(206) = EXP(+(-7.23838D+05)*TM3+(-2.77495D+04)*TM2                      
     &          +(-5.83020D+03)*TM1+( 3.04837D+01)+( 5.16215D-04)*TP1           
     &          +(-1.10029D-07)*TP2+( 9.87866D-12)*TP3)                         
      RF(207) = 5.0000D+13                                                      
      RB(207) = EXP(+( 2.81353D+06)*TM3+(-7.21656D+04)*TM2                      
     &          +(-5.57802D+04)*TM1+( 3.25521D+01)+( 8.58276D-04)*TP1           
     &          +(-1.82773D-07)*TP2+( 1.65515D-11)*TP3)                         
      RF(208) = 2.0000D+13                                                      
      RB(208) = EXP(+(-7.95066D+06)*TM3+( 5.81745D+04)*TM2                      
     &          +(-6.39100D+04)*TM1+( 3.37276D+01)+( 6.50408D-04)*TP1           
     &          +(-1.71005D-07)*TP2+( 1.72129D-11)*TP3)                         
      RF(209) = 2.5000D+13                                                      
      RB(209) = EXP(+(-3.28303D+07)*TM3+( 3.37854D+05)*TM2                      
     &          +(-5.78289D+04)*TM1+( 3.71230D+01)+( 4.99577D-05)*TP1           
     &          +(-1.12104D-07)*TP2+( 1.36151D-11)*TP3)                         
      RF(210) = 8.9500D+19 * EXP( -1.320*ALOGT -  0.74008*RTR )                 
      RB(210) = EXP(+(-6.78059D+05)*TM3+( 1.61997D+04)*TM2                      
     &          +(-2.42505D+04)*TM1+( 3.91620D+01)+(-7.81019D-04)*TP1           
     &          +(-4.07714D-09)*TP2+( 5.44731D-12)*TP3)                         
      RF(211) = 2.5000D+13                                                      
      RB(211) = EXP(+( 6.22282D+06)*TM3+(-9.12968D+04)*TM2                      
     &          +(-2.66121D+04)*TM1+( 2.86193D+01)+( 2.76536D-04)*TP1           
     &          +( 2.33139D-08)*TP2+(-4.51236D-12)*TP3)                         
      RF(212) = 4.5000D+11 * EXP(  0.720*ALOGT -  0.66007*RTR )                 
      RB(212) = EXP(+( 1.55806D+06)*TM3+(-3.14295D+04)*TM2                      
     &          +(-2.83427D+04)*TM1+( 3.02472D+01)+( 8.45949D-04)*TP1           
     &          +(-6.90876D-08)*TP2+( 2.63570D-12)*TP3)                         
      RF(213) = 1.3000D+07 * EXP(  1.900*ALOGT +  0.95010*RTR )                 
      RB(213) = EXP(+(-2.07156D+07)*TM3+( 2.36532D+05)*TM2                      
     &          +(-3.64866D+04)*TM1+( 3.00194D+01)+( 1.48984D-03)*TP1           
     &          +(-2.05248D-07)*TP2+( 1.63082D-11)*TP3)                         
      RF(214) = 1.0000D+13 * EXP(              - 13.00141*RTR )                 
      RB(214) = EXP(+( 7.50241D+06)*TM3+(-7.42266D+04)*TM2                      
     &          +(-6.58743D+03)*TM1+( 2.80025D+01)+(-5.26844D-05)*TP1           
     &          +( 7.80162D-08)*TP2+(-9.47738D-12)*TP3)                         
      RF(215) = 7.7000D+13                                                      
      RB(215) = EXP(+(-3.12079D+06)*TM3+( 3.10132D+04)*TM2                      
     &          +(-3.92807D+04)*TM1+( 3.36363D+01)+( 1.26958D-06)*TP1           
     &          +(-1.06463D-08)*TP2+( 4.21622D-12)*TP3)                         
      RF(216) = 4.0000D+13                                                      
      RB(216) = EXP(+( 8.56992D+06)*TM3+(-4.37460D+04)*TM2                      
     &          +(-1.55299D+04)*TM1+( 3.70081D+01)+(-1.56081D-03)*TP1           
     &          +( 3.36044D-07)*TP2+(-2.94961D-11)*TP3)                         
      RF(217) = 8.0000D+12 * EXP(              -  7.46081*RTR )                 
      RB(217) = EXP(+( 1.08105D+07)*TM3+(-1.00820D+05)*TM2                      
     &          +(-6.19404D+03)*TM1+( 3.04213D+01)+(-6.09751D-04)*TP1           
     &          +( 1.16547D-07)*TP2+(-7.80214D-12)*TP3)                         
      RF(218) = 6.1400D+12 * EXP(              +  0.44005*RTR )                 
      RB(218) = EXP(+( 1.13670D+07)*TM3+(-9.89721D+04)*TM2                      
     &          +(-6.50067D+03)*TM1+( 3.13665D+01)+(-8.68668D-04)*TP1           
     &          +( 1.65263D-07)*TP2+(-1.18494D-11)*TP3)                         
      RF(219) = 2.1000D+13 * EXP(              -  4.71051*RTR )                 
      RB(219) = EXP(+( 4.63345D+04)*TM3+( 2.95202D+04)*TM2                      
     &          +(-1.29399D+04)*TM1+( 3.34781D+01)+(-8.17619D-04)*TP1           
     &          +( 1.28314D-07)*TP2+(-7.14070D-12)*TP3)                         
      RF(220) = 2.3500D+13                                                      
      RB(220) = EXP(+(-1.69056D+07)*TM3+( 1.38430D+05)*TM2                      
     &          +(-4.84326D+04)*TM1+( 2.87850D+01)+( 9.81242D-04)*TP1           
     &          +(-2.05369D-07)*TP2+( 2.00074D-11)*TP3)                         
      RF(221) = 5.4000D+13                                                      
      RB(221) = EXP(+(-1.41917D+07)*TM3+( 9.64699D+04)*TM2                      
     &          +(-1.25089D+04)*TM1+( 2.70838D+01)+( 1.45217D-03)*TP1           
     &          +(-3.24045D-07)*TP2+( 3.14014D-11)*TP3)                         
      RF(222) = 2.5000D+12                                                      
      RB(222) = EXP(+(-9.57534D+06)*TM3+( 5.95774D+04)*TM2                      
     &          +( 2.97964D+03)*TM1+( 2.64563D+01)+( 5.32909D-04)*TP1           
     &          +(-5.91270D-08)*TP2+( 4.51771D-12)*TP3)                         
      RF(223) = 2.0000D+13                                                      
      RB(223) = EXP(+(-1.82271D+07)*TM3+( 1.39829D+05)*TM2                      
     &          +(-8.62452D+04)*TM1+( 2.98361D+01)+( 1.16370D-03)*TP1           
     &          +(-2.53327D-07)*TP2+( 2.51601D-11)*TP3)                         
      RF(224) = 2.0000D+12 * EXP(              - 20.00217*RTR )                 
      RB(224) = EXP(+(-7.48196D+06)*TM3+( 7.05778D+04)*TM2                      
     &          +(-6.24750D+04)*TM1+( 2.99896D+01)+( 2.90958D-04)*TP1           
     &          +(-8.74075D-08)*TP2+( 1.06863D-11)*TP3)                         
      RF(225) = 8.8000D+16 * EXP( -0.500*ALOGT - 48.00521*RTR )                 
      RB(225) = EXP(+( 5.16394D+05)*TM3+(-6.24070D+04)*TM2                      
     &          +( 3.85368D+03)*TM1+( 3.14553D+01)+( 7.52826D-04)*TP1           
     &          +(-1.37768D-07)*TP2+( 1.27095D-11)*TP3)                         
      RF(226) = 2.8500D+17 * EXP( -1.520*ALOGT -  0.74008*RTR )                 
      RB(226) = EXP(+( 1.13862D+07)*TM3+(-1.43663D+05)*TM2                      
     &          +(-2.97160D+04)*TM1+( 3.27961D+01)+(-8.69877D-04)*TP1           
     &          +( 1.18884D-07)*TP2+(-8.95848D-12)*TP3)                         
      RF(227) = 5.7000D+18 * EXP( -2.000*ALOGT -  0.80009*RTR )                 
      RB(227) = EXP(+( 1.31217D+07)*TM3+(-1.69723D+05)*TM2                      
     &          +(-7.33537D+04)*TM1+( 3.40246D+01)+(-1.08155D-03)*TP1           
     &          +( 1.44820D-07)*TP2+(-1.01553D-11)*TP3)                         
      RF(228) = 1.0400D+29 * EXP( -3.300*ALOGT -126.61374*RTR )                 
      RB(228) = EXP(+( 3.83330D+08)*TM3+(-3.76427D+07)*TM2                      
     &          +( 1.63247D+05)*TM1+(-2.17848D+02)+( 1.85741D-01)*TP1           
     &          +(-6.22578D-05)*TP2+( 7.73448D-09)*TP3)                         
      RF(229) = 1.1070D+04 * EXP(  2.640*ALOGT -  4.98054*RTR )                 
      RB(229) = EXP(+(-1.95843D+07)*TM3+( 2.58737D+05)*TM2                      
     &          +(-8.41541D+03)*TM1+( 2.97094D+01)+( 1.11275D-03)*TP1           
     &          +(-1.21087D-07)*TP2+( 7.54492D-12)*TP3)                         
      RF(230) = 2.7670D+03 * EXP(  2.640*ALOGT -  4.98054*RTR )                 
      RB(230) = EXP(+(-3.37760D+07)*TM3+( 3.55207D+05)*TM2                      
     &          +(-2.09244D+04)*TM1+( 2.37867D+01)+( 2.56492D-03)*TP1           
     &          +(-4.45131D-07)*TP2+( 3.89463D-11)*TP3)                         
      RF(231) = 2.1340D+09 * EXP(  1.580*ALOGT - 26.60289*RTR )                 
      RB(231) = EXP(+(-1.78152D+07)*TM3+( 1.78857D+05)*TM2                      
     &          +(-3.02577D+03)*TM1+( 2.88849D+01)+( 1.90842D-03)*TP1           
     &          +(-3.24247D-07)*TP2+( 2.53532D-11)*TP3)                         
      RF(232) = 1.1000D+06 * EXP(  2.030*ALOGT - 13.37145*RTR )                 
      RB(232) = EXP(+(-2.45959D+07)*TM3+( 3.00605D+05)*TM2                      
     &          +(-4.66439D+03)*TM1+( 3.22884D+01)+( 3.77984D-04)*TP1           
     &          +( 1.92403D-08)*TP2+(-4.12534D-12)*TP3)                         
      RF(233) = 4.4000D+03 * EXP(  2.260*ALOGT -  6.40069*RTR )                 
      RB(233) = EXP(+(-2.85166D+07)*TM3+( 3.52723D+05)*TM2                      
     &          +(-1.42430D+04)*TM1+( 2.88921D+01)+( 1.40449D-04)*TP1           
     &          +( 7.57252D-08)*TP2+(-1.02337D-11)*TP3)                         
      RF(234) = 1.6000D+02 * EXP(  2.560*ALOGT -  9.00098*RTR )                 
      RB(234) = EXP(+(-4.30584D+07)*TM3+( 4.67436D+05)*TM2                      
     &          +(-1.79939D+04)*TM1+( 2.22881D+01)+( 2.11968D-03)*TP1           
     &          +(-3.74898D-07)*TP2+( 3.31966D-11)*TP3)                         
      RF(235) = 1.4000D+26 * EXP( -3.400*ALOGT -  1.90021*RTR )                 
      RB(235) = EXP(+( 1.12448D+07)*TM3+(-1.19020D+05)*TM2                      
     &          +(-1.21586D+04)*TM1+( 3.85402D+01)+(-2.99971D-03)*TP1           
     &          +( 5.23347D-07)*TP2+(-4.69309D-11)*TP3)                         
      RF(236) = 6.0000D+13 * EXP(              -  0.40004*RTR )                 
      RB(236) = EXP(+( 2.86881D+06)*TM3+(-1.00009D+05)*TM2                      
     &          +(-3.89565D+04)*TM1+( 2.91973D+01)+( 1.21687D-03)*TP1           
     &          +(-2.77008D-07)*TP2+( 2.89732D-11)*TP3)                         
      RF(237) = 6.3000D+13 * EXP(              - 46.02499*RTR )                 
      RB(237) = EXP(+( 2.51630D+06)*TM3+(-1.64365D+04)*TM2                      
     &          +( 2.78415D+02)*TM1+( 3.10434D+01)+(-8.08340D-05)*TP1           
     &          +( 2.81445D-08)*TP2+(-5.41201D-12)*TP3)                         
      RF(238) = 2.8570D+08 * EXP(  1.100*ALOGT - 20.40221*RTR )                 
      RB(238) = EXP(+(-2.54535D+06)*TM3+( 8.82789D+04)*TM2                      
     &          +(-9.84829D+03)*TM1+( 2.99309D+01)+( 3.80111D-05)*TP1           
     &          +( 3.24882D-08)*TP2+(-4.36883D-12)*TP3)                         
      RF(240) = 1.0000D+13 * EXP(              - 74.00803*RTR )                 
      RB(240) = EXP(+( 9.21797D+06)*TM3+(-7.69439D+04)*TM2                      
     &          +(-2.50798D+04)*TM1+( 3.02711D+01)+(-3.32956D-04)*TP1           
     &          +( 5.62979D-08)*TP2+(-6.35304D-12)*TP3)                         
      RF(241) = 1.0000D+11 * EXP(              - 65.00705*RTR )                 
      RB(241) = EXP(+( 1.81274D+07)*TM3+(-1.65062D+05)*TM2                      
     &          +(-2.47665D+04)*TM1+( 2.46659D+01)+(-2.59327D-04)*TP1           
     &          +( 4.30146D-08)*TP2+(-5.44344D-12)*TP3)                         
      RF(242) = 1.9000D+13                                                      
      RB(242) = EXP(+( 1.19472D+06)*TM3+(-1.50368D+04)*TM2                      
     &          +(-1.43737D+04)*TM1+( 3.10571D+01)+( 1.01627D-04)*TP1           
     &          +(-1.98132D-08)*TP2+(-2.59253D-13)*TP3)                         
      RF(243) = 2.9000D+13                                                      
      RB(243) = EXP(+(-1.92608D+06)*TM3+( 1.59764D+04)*TM2                      
     &          +(-5.36544D+04)*TM1+( 3.31414D+01)+( 1.02896D-04)*TP1           
     &          +(-3.04595D-08)*TP2+( 3.95696D-12)*TP3)                         
      RF(244) = 5.0000D+13                                                      
      RB(244) = EXP(+( 6.86221D+06)*TM3+(-3.86124D+04)*TM2                      
     &          +(-3.66258D+04)*TM1+( 3.56346D+01)+(-5.73542D-04)*TP1           
     &          +( 1.22429D-07)*TP2+(-1.13450D-11)*TP3)                         
      RF(245) = 2.0000D+13                                                      
      RB(245) = EXP(+( 1.30278D+07)*TM3+(-8.77735D+04)*TM2                      
     &          +(-4.06907D+04)*TM1+( 3.69209D+01)+(-1.36642D-03)*TP1           
     &          +( 3.32300D-07)*TP2+(-3.29096D-11)*TP3)                         
      RF(246) = 3.0000D+13                                                      
      RB(246) = EXP(+(-5.94961D+06)*TM3+( 9.16117D+04)*TM2                      
     &          +(-2.13782D+04)*TM1+( 3.38817D+01)+(-6.16922D-04)*TP1           
     &          +( 1.60229D-07)*TP2+(-1.53208D-11)*TP3)                         
      RF(247) = 3.1000D+17 * EXP( -1.380*ALOGT -  1.27014*RTR )                 
      RB(247) = EXP(+( 1.73846D+07)*TM3+(-1.69059D+05)*TM2                      
     &          +(-4.60107D+04)*TM1+( 3.77505D+01)+(-2.52761D-03)*TP1           
     &          +( 5.17740D-07)*TP2+(-4.92588D-11)*TP3)                         
      RF(248) = 2.9000D+14 * EXP( -0.690*ALOGT -  0.76008*RTR )                 
      RB(248) = EXP(+( 1.71275D+07)*TM3+(-1.77728D+05)*TM2                      
     &          +(-3.67927D+04)*TM1+( 3.06047D+01)+(-5.38653D-04)*TP1           
     &          +( 7.21946D-08)*TP2+(-6.49751D-12)*TP3)                         
      RF(249) = 3.8000D+13 * EXP( -0.360*ALOGT -  0.58006*RTR )                 
      RB(249) = EXP(+( 1.65246D+07)*TM3+(-1.27679D+05)*TM2                      
     &          +(-1.13933D+04)*TM1+( 3.53967D+01)+(-1.91602D-03)*TP1           
     &          +( 4.14076D-07)*TP2+(-4.00240D-11)*TP3)                         
      RF(250) = 3.1000D+17 * EXP( -1.380*ALOGT -  1.27014*RTR )                 
      RB(250) = EXP(+( 2.62940D+07)*TM3+(-2.57177D+05)*TM2                      
     &          +(-5.02269D+04)*TM1+( 3.67505D+01)+(-2.45398D-03)*TP1           
     &          +( 5.04456D-07)*TP2+(-4.83492D-11)*TP3)                         
      RF(251) = 2.9000D+14 * EXP( -0.690*ALOGT -  0.76008*RTR )                 
      RB(251) = EXP(+( 2.60369D+07)*TM3+(-2.65846D+05)*TM2                      
     &          +(-4.10089D+04)*TM1+( 2.96047D+01)+(-4.65024D-04)*TP1           
     &          +( 5.89113D-08)*TP2+(-5.58791D-12)*TP3)                         
      RF(252) = 3.8000D+13 * EXP( -0.360*ALOGT -  0.58006*RTR )                 
      RB(252) = EXP(+( 2.54341D+07)*TM3+(-2.15797D+05)*TM2                      
     &          +(-1.56095D+04)*TM1+( 3.43967D+01)+(-1.84239D-03)*TP1           
     &          +( 4.00793D-07)*TP2+(-3.91144D-11)*TP3)                         
      RF(253) = 9.6000D+13 * EXP(              - 28.80312*RTR )                 
      RB(253) = EXP(+( 1.51255D+07)*TM3+(-1.53788D+05)*TM2                      
     &          +(-5.54232D+04)*TM1+( 3.30314D+01)+( 4.64672D-04)*TP1           
     &          +(-1.09840D-07)*TP2+( 8.97707D-12)*TP3)                         
      RF(254) = 1.0000D+12 * EXP(              - 21.75236*RTR )                 
      RB(254) = EXP(+( 8.94391D+06)*TM3+(-6.13589D+04)*TM2                      
     &          +(-5.01502D+03)*TM1+( 2.73733D+01)+(-3.70975D-04)*TP1           
     &          +( 1.23889D-07)*TP2+(-1.58246D-11)*TP3)                         
      RF(255) = 2.2000D+13                                                      
      RB(255) = EXP(+(-1.64526D+07)*TM3+( 1.05288D+05)*TM2                      
     &          +(-7.29814D+04)*TM1+( 2.99297D+01)+( 1.35351D-03)*TP1           
     &          +(-2.39438D-07)*TP2+( 2.16151D-11)*TP3)                         
      RF(256) = 2.0000D+12                                                      
      RB(256) = EXP(+(-1.69361D+07)*TM3+( 1.57143D+05)*TM2                      
     &          +(-5.87653D+04)*TM1+( 2.82317D+01)+( 8.50191D-04)*TP1           
     &          +(-2.00902D-07)*TP2+( 1.98042D-11)*TP3)                         
      RF(257) = 1.2000D+13                                                      
      RB(257) = EXP(+(-2.09422D+07)*TM3+( 1.54688D+05)*TM2                      
     &          +(-2.11104D+04)*TM1+( 2.61381D+01)+( 1.23308D-03)*TP1           
     &          +(-2.35599D-07)*TP2+( 2.31381D-11)*TP3)                         
      RF(258) = 1.2000D+13                                                      
      RB(258) = EXP(+(-2.37393D+07)*TM3+( 2.09915D+05)*TM2                      
     &          +(-2.99182D+04)*TM1+( 2.99056D+01)+( 5.40939D-04)*TP1           
     &          +(-6.48178D-08)*TP2+( 5.49144D-12)*TP3)                         
      RF(259) = 1.0000D+14                                                      
      RB(259) = EXP(+(-2.34402D+07)*TM3+( 1.92127D+05)*TM2                      
     &          +(-3.50039D+04)*TM1+( 2.92731D+01)+( 1.65407D-03)*TP1           
     &          +(-3.75876D-07)*TP2+( 3.75513D-11)*TP3)                         
      RF(260) = 9.8000D+07 * EXP(  1.410*ALOGT -  8.50092*RTR )                 
      RB(260) = EXP(+(-8.67940D+06)*TM3+( 1.09986D+05)*TM2                      
     &          +(-2.51781D+04)*TM1+( 2.82937D+01)+( 1.78554D-03)*TP1           
     &          +(-3.61237D-07)*TP2+( 3.35694D-11)*TP3)                         
      RF(261) = 1.5000D+08 * EXP(  1.570*ALOGT - 44.00477*RTR )                 
      RB(261) = EXP(+(-2.58031D+07)*TM3+( 2.74528D+05)*TM2                      
     &          +(-3.96828D+04)*TM1+( 2.71500D+01)+( 2.53599D-03)*TP1           
     &          +(-5.74676D-07)*TP2+( 5.54197D-11)*TP3)                         
      RF(262) = 2.2000D+06 * EXP(  2.110*ALOGT - 11.40124*RTR )                 
      RB(262) = EXP(+(-7.94181D+06)*TM3+( 1.07781D+05)*TM2                      
     &          +(-1.83285D+03)*TM1+( 2.64218D+01)+( 2.22107D-03)*TP1           
     &          +(-4.13689D-07)*TP2+( 3.68542D-11)*TP3)                         
      RF(263) = 2.2500D+07 * EXP(  1.700*ALOGT -  3.80041*RTR )                 
      RB(263) = EXP(+(-2.81972D+07)*TM3+( 2.77993D+05)*TM2                      
     &          +(-5.33276D+03)*TM1+( 2.32881D+01)+( 2.98980D-03)*TP1           
     &          +(-6.26131D-07)*TP2+( 5.88671D-11)*TP3)                         
      RF(264) = 1.0500D+05 * EXP(  2.500*ALOGT - 13.30144*RTR )                 
      RB(264) = EXP(+(-9.38783D+06)*TM3+( 1.29161D+05)*TM2                      
     &          +(-3.95702D+03)*TM1+( 2.67505D+01)+( 2.55228D-03)*TP1           
     &          +(-4.64721D-07)*TP2+( 4.03635D-11)*TP3)                         
      RF(265) = 4.6500D+12 * EXP(              -  6.85074*RTR )                 
      RB(265) = EXP(+( 4.23242D+06)*TM3+(-3.20691D+04)*TM2                      
     &          +(-7.09434D+03)*TM1+( 2.92185D+01)+( 5.39835D-04)*TP1           
     &          +(-1.39547D-07)*TP2+( 1.34592D-11)*TP3)                         
      RF(266) = 1.5500D+12 * EXP(              -  6.85074*RTR )                 
      RB(266) = EXP(+(-4.98927D+06)*TM3+( 6.70994D+04)*TM2                      
     &          +(-1.84649D+04)*TM1+( 3.01479D+01)+( 3.80264D-04)*TP1           
     &          +(-1.24272D-07)*TP2+( 1.31547D-11)*TP3)                         
      RF(267) = 1.1800D+16 * EXP(              - 84.72919*RTR )                 
      RB(267) = EXP(+(-9.61549D+08)*TM3+( 6.61632D+06)*TM2                      
     &          +(-1.59127D+04)*TM1+( 5.14953D+01)+(-1.20264D-02)*TP1           
     &          +( 3.85132D-06)*TP2+(-4.57491D-10)*TP3)                         
      RF(268) = 2.1000D+15 * EXP( -0.690*ALOGT -  2.85031*RTR )                 
      RB(268) = EXP(+(-2.35878D+06)*TM3+(-2.89252D+03)*TM2                      
     &          +(-3.59350D+04)*TM1+( 3.09021D+01)+(-3.73385D-04)*TP1           
     &          +( 6.22934D-08)*TP2+(-5.59618D-12)*TP3)                         
      RF(269) = 2.7000D+11 * EXP(  0.180*ALOGT -  2.12023*RTR )                 
      RB(269) = EXP(+(-4.37154D+06)*TM3+( 9.43108D+03)*TM2                      
     &          +(-2.67320D+04)*TM1+( 2.30123D+01)+( 1.74550D-03)*TP1           
     &          +(-4.05817D-07)*TP2+( 3.91499D-11)*TP3)                         
      RF(270) = 1.7000D+14 * EXP( -0.750*ALOGT -  2.89031*RTR )                 
      RB(270) = EXP(+(-1.33893D+07)*TM3+( 6.98348D+04)*TM2                      
     &          +(-3.81460D+04)*TM1+( 2.26174D+01)+( 1.34599D-03)*TP1           
     &          +(-3.43200D-07)*TP2+( 3.38646D-11)*TP3)                         
      RF(271) = 2.0000D+07 * EXP(  2.000*ALOGT -  2.00022*RTR )                 
      RB(271) = EXP(+(-2.11848D+07)*TM3+( 2.58550D+05)*TM2                      
     &          +(-1.53294D+04)*TM1+( 3.11351D+01)+( 1.04011D-03)*TP1           
     &          +(-1.65407D-07)*TP2+( 1.34082D-11)*TP3)                         
      RF(272) = 2.3500D+13                                                      
      RB(272) = EXP(+(-9.24567D+06)*TM3+( 1.02017D+05)*TM2                      
     &          +(-2.54359D+04)*TM1+( 3.34607D+01)+(-2.79713D-04)*TP1           
     &          +( 6.24081D-08)*TP2+(-5.47609D-12)*TP3)                         
      RF(273) = 6.1000D+14 * EXP( -0.310*ALOGT -  0.29003*RTR )                 
      RB(273) = EXP(+( 6.75274D+06)*TM3+(-3.38433D+04)*TM2                      
     &          +(-1.86800D+04)*TM1+( 3.36741D+01)+(-1.17557D-03)*TP1           
     &          +( 3.04073D-07)*TP2+(-3.29477D-11)*TP3)                         
      RF(274) = 3.7000D+12 * EXP(  0.150*ALOGT +  0.09001*RTR )                 
      RB(274) = EXP(+( 1.92117D+07)*TM3+(-2.02963D+05)*TM2                      
     &          +(-5.75404D+04)*TM1+( 3.07411D+01)+( 1.99982D-04)*TP1           
     &          +( 9.09739D-10)*TP2+(-3.73538D-12)*TP3)                         
      RF(275) = 5.4000D+05 * EXP(  2.400*ALOGT -  9.91608*RTR )                 
      RB(275) = EXP(+(-1.23125D+07)*TM3+( 1.23351D+05)*TM2                      
     &          +(-3.54788D+03)*TM1+( 2.63302D+01)+( 2.26829D-03)*TP1           
     &          +(-3.75962D-07)*TP2+( 3.05199D-11)*TP3)                         
      RF(276) = 5.0000D+07 * EXP(  1.600*ALOGT -  0.95510*RTR )                 
      RB(276) = EXP(+(-1.52737D+07)*TM3+( 1.60388D+05)*TM2                      
     &          +(-6.60954D+03)*TM1+( 2.74366D+01)+( 1.48296D-03)*TP1           
     &          +(-2.63905D-07)*TP2+( 2.23603D-11)*TP3)                         
      RF(277) = 9.4000D+06 * EXP(  1.940*ALOGT -  6.46070*RTR )                 
      RB(277) = EXP(+(-1.01837D+07)*TM3+( 9.38068D+04)*TM2                      
     &          +(-5.92221D+02)*TM1+( 2.53337D+01)+( 1.88656D-03)*TP1           
     &          +(-3.16155D-07)*TP2+( 2.62387D-11)*TP3)                         
      TOLD=T
      END IF
      RFH49 = 2.5000D+16 * EXP( -0.800*ALOGT                )                   
      RFL49 = 3.2000D+27 * EXP( -3.140*ALOGT -  1.23013*RTR )                   
      PR = RFL49*XM(49)/RFH49                                                   
      PRLOG = DLOG10(MAX(PR,SMALL))
      F4 =     6.80000D-01
      F5 =     7.80000D+01
      F6 =     1.99500D+03
      F7 =     5.59000D+03
      FC = (1.0D0-F4)*EXP(-T/F5)+F4*EXP(-T/F6)+EXP(-F7/T)
      FCLOG = DLOG10(MAX(FC,SMALL))
      CPRLOG = PRLOG - (0.4D0 + 0.67D0*FCLOG)
      X = CPRLOG/(0.75D0-1.27D0*FCLOG-0.14D0*CPRLOG)
      F = 10.0D0**( FCLOG/(1.0D0+X*X) )
      PCOR = PR*F/(1.0D0+PR)
      RF(49) = RFH49*PCOR                                                       
      RB(49) = EXP(+(-1.26615D+07)*TM3+( 1.48317D+05)*TM2                       
     &         +(-5.58420D+04)*TM1+( 3.63302D+01)+(-7.92566D-04)*TP1            
     &         +( 5.92099D-08)*TP2+(-1.32724D-12)*TP3)                          
      RB(49) = RB(49)*PCOR                                                      
      RFH51 = 1.2700D+16 * EXP( -0.630*ALOGT -  0.38304*RTR )                   
      RFL51 = 2.4770D+33 * EXP( -4.760*ALOGT -  2.44026*RTR )                   
      PR = RFL51*XM(51)/RFH51                                                   
      PRLOG = DLOG10(MAX(PR,SMALL))
      F4 =     7.83000D-01
      F5 =     7.40000D+01
      F6 =     2.94100D+03
      F7 =     6.96400D+03
      FC = (1.0D0-F4)*EXP(-T/F5)+F4*EXP(-T/F6)+EXP(-F7/T)
      FCLOG = DLOG10(MAX(FC,SMALL))
      CPRLOG = PRLOG - (0.4D0 + 0.67D0*FCLOG)
      X = CPRLOG/(0.75D0-1.27D0*FCLOG-0.14D0*CPRLOG)
      F = 10.0D0**( FCLOG/(1.0D0+X*X) )
      PCOR = PR*F/(1.0D0+PR)
      RF(51) = RFH51*PCOR                                                       
      RB(51) = EXP(+(-3.44713D+07)*TM3+( 3.91291D+05)*TM2                       
     &         +(-5.41091D+04)*TM1+( 3.87736D+01)+(-7.65043D-04)*TP1            
     &         +( 1.26478D-09)*TP2+( 4.81589D-12)*TP3)                          
      RB(51) = RB(51)*PCOR                                                      
      RFH53 = 1.0900D+12 * EXP(  0.480*ALOGT +  0.26003*RTR )                   
      RFL53 = 1.3500D+24 * EXP( -2.570*ALOGT -  1.42515*RTR )                   
      PR = RFL53*XM(53)/RFH53                                                   
      PRLOG = DLOG10(MAX(PR,SMALL))
      F4 =     7.82400D-01
      F5 =     2.71000D+02
      F6 =     2.75500D+03
      F7 =     6.57000D+03
      FC = (1.0D0-F4)*EXP(-T/F5)+F4*EXP(-T/F6)+EXP(-F7/T)
      FCLOG = DLOG10(MAX(FC,SMALL))
      CPRLOG = PRLOG - (0.4D0 + 0.67D0*FCLOG)
      X = CPRLOG/(0.75D0-1.27D0*FCLOG-0.14D0*CPRLOG)
      F = 10.0D0**( FCLOG/(1.0D0+X*X) )
      PCOR = PR*F/(1.0D0+PR)
      RF(53) = RFH53*PCOR                                                       
      RB(53) = EXP(+(-2.96617D+07)*TM3+( 3.57134D+05)*TM2                       
     &         +(-4.54461D+04)*TM1+( 3.60768D+01)+( 1.04094D-04)*TP1            
     &         +(-9.93140D-08)*TP2+( 1.31708D-11)*TP3)                          
      RB(53) = RB(53)*PCOR                                                      
      RFH55 = 5.4000D+11 * EXP(  0.454*ALOGT -  3.60039*RTR )                   
      RFL55 = 1.2700D+32 * EXP( -4.820*ALOGT -  6.53071*RTR )                   
      PR = RFL55*XM(55)/RFH55                                                   
      PRLOG = DLOG10(MAX(PR,SMALL))
      F4 =     7.18700D-01
      F5 =     1.03000D+02
      F6 =     1.29100D+03
      F7 =     4.16000D+03
      FC = (1.0D0-F4)*EXP(-T/F5)+F4*EXP(-T/F6)+EXP(-F7/T)
      FCLOG = DLOG10(MAX(FC,SMALL))
      CPRLOG = PRLOG - (0.4D0 + 0.67D0*FCLOG)
      X = CPRLOG/(0.75D0-1.27D0*FCLOG-0.14D0*CPRLOG)
      F = 10.0D0**( FCLOG/(1.0D0+X*X) )
      PCOR = PR*F/(1.0D0+PR)
      RF(55) = RFH55*PCOR                                                       
      RB(55) = EXP(+( 6.31058D+06)*TM3+(-4.77943D+04)*TM2                       
     &         +(-1.64365D+04)*TM1+( 2.94949D+01)+( 3.14782D-04)*TP1            
     &         +(-1.01093D-07)*TP2+( 1.04005D-11)*TP3)                          
      RB(55) = RB(55)*PCOR                                                      
      RFH56 = 5.4000D+11 * EXP(  0.454*ALOGT -  2.60028*RTR )                   
      RFL56 = 2.2000D+30 * EXP( -4.800*ALOGT -  5.56060*RTR )                   
      PR = RFL56*XM(56)/RFH56                                                   
      PRLOG = DLOG10(MAX(PR,SMALL))
      F4 =     7.58000D-01
      F5 =     9.40000D+01
      F6 =     1.55500D+03
      F7 =     4.20000D+03
      FC = (1.0D0-F4)*EXP(-T/F5)+F4*EXP(-T/F6)+EXP(-F7/T)
      FCLOG = DLOG10(MAX(FC,SMALL))
      CPRLOG = PRLOG - (0.4D0 + 0.67D0*FCLOG)
      X = CPRLOG/(0.75D0-1.27D0*FCLOG-0.14D0*CPRLOG)
      F = 10.0D0**( FCLOG/(1.0D0+X*X) )
      PCOR = PR*F/(1.0D0+PR)
      RF(56) = RFH56*PCOR                                                       
      RB(56) = EXP(+(-1.72463D+07)*TM3+( 2.38003D+05)*TM2                       
     &         +(-1.33724D+04)*TM1+( 3.30241D+01)+(-7.76132D-05)*TP1            
     &         +(-2.92251D-08)*TP2+( 2.95920D-12)*TP3)                          
      RB(56) = RB(56)*PCOR                                                      
      RFH58 = 1.8000D+13                                                        
      RFL58 = 3.0000D+31 * EXP( -4.800*ALOGT -  3.30036*RTR )                   
      PR = RFL58*XM(58)/RFH58                                                   
      PRLOG = DLOG10(MAX(PR,SMALL))
      F4 =     7.67900D-01
      F5 =     3.38000D+02
      F6 =     1.81200D+03
      F7 =     5.08100D+03
      FC = (1.0D0-F4)*EXP(-T/F5)+F4*EXP(-T/F6)+EXP(-F7/T)
      FCLOG = DLOG10(MAX(FC,SMALL))
      CPRLOG = PRLOG - (0.4D0 + 0.67D0*FCLOG)
      X = CPRLOG/(0.75D0-1.27D0*FCLOG-0.14D0*CPRLOG)
      F = 10.0D0**( FCLOG/(1.0D0+X*X) )
      PCOR = PR*F/(1.0D0+PR)
      RF(58) = RFH58*PCOR                                                       
      RB(58) = EXP(+(-3.60222D+07)*TM3+( 4.34426D+05)*TM2                       
     &         +(-5.00604D+04)*TM1+( 3.60552D+01)+(-2.15830D-04)*TP1            
     &         +(-4.81374D-08)*TP2+( 8.80885D-12)*TP3)                          
      RB(58) = RB(58)*PCOR                                                      
      RFH62 = 5.0000D+13                                                        
      RFL62 = 8.6000D+28 * EXP( -4.000*ALOGT -  3.02533*RTR )                   
      PR = RFL62*XM(62)/RFH62                                                   
      PRLOG = DLOG10(MAX(PR,SMALL))
      F4 =     8.90200D-01
      F5 =     1.44000D+02
      F6 =     2.83800D+03
      F7 =     4.55690D+04
      FC = (1.0D0-F4)*EXP(-T/F5)+F4*EXP(-T/F6)+EXP(-F7/T)
      FCLOG = DLOG10(MAX(FC,SMALL))
      CPRLOG = PRLOG - (0.4D0 + 0.67D0*FCLOG)
      X = CPRLOG/(0.75D0-1.27D0*FCLOG-0.14D0*CPRLOG)
      F = 10.0D0**( FCLOG/(1.0D0+X*X) )
      PCOR = PR*F/(1.0D0+PR)
      RF(62) = RFH62*PCOR                                                       
      RB(62) = EXP(+(-1.24653D+07)*TM3+( 1.48628D+05)*TM2                       
     &         +(-5.26213D+04)*TM1+( 3.35477D+01)+( 1.76565D-04)*TP1            
     &         +(-1.20006D-07)*TP2+( 1.62501D-11)*TP3)                          
      RB(62) = RB(62)*PCOR                                                      
      RFH69 = 1.0000D+17 * EXP( -1.000*ALOGT                )                   
      RFL69 = 3.7500D+33 * EXP( -4.800*ALOGT -  1.90021*RTR )                   
      PR = RFL69*XM(69)/RFH69                                                   
      PRLOG = DLOG10(MAX(PR,SMALL))
      F4 =     6.46400D-01
      F5 =     1.32000D+02
      F6 =     1.31500D+03
      F7 =     5.56600D+03
      FC = (1.0D0-F4)*EXP(-T/F5)+F4*EXP(-T/F6)+EXP(-F7/T)
      FCLOG = DLOG10(MAX(FC,SMALL))
      CPRLOG = PRLOG - (0.4D0 + 0.67D0*FCLOG)
      X = CPRLOG/(0.75D0-1.27D0*FCLOG-0.14D0*CPRLOG)
      F = 10.0D0**( FCLOG/(1.0D0+X*X) )
      PCOR = PR*F/(1.0D0+PR)
      RF(69) = RFH69*PCOR                                                       
      RB(69) = EXP(+( 7.67516D+06)*TM3+(-2.41211D+04)*TM2                       
     &         +(-6.64460D+04)*TM1+( 3.74795D+01)+(-1.02498D-03)*TP1            
     &         +( 1.57486D-07)*TP2+(-1.48357D-11)*TP3)                          
      RB(69) = RB(69)*PCOR                                                      
      RFH70 = 5.6000D+12 * EXP(              -  2.40026*RTR )                   
      RFL70 = 3.8000D+40 * EXP( -7.270*ALOGT -  7.22078*RTR )                   
      PR = RFL70*XM(70)/RFH70                                                   
      PRLOG = DLOG10(MAX(PR,SMALL))
      F4 =     7.50700D-01
      F5 =     9.85000D+01
      F6 =     1.30200D+03
      F7 =     4.16700D+03
      FC = (1.0D0-F4)*EXP(-T/F5)+F4*EXP(-T/F6)+EXP(-F7/T)
      FCLOG = DLOG10(MAX(FC,SMALL))
      CPRLOG = PRLOG - (0.4D0 + 0.67D0*FCLOG)
      X = CPRLOG/(0.75D0-1.27D0*FCLOG-0.14D0*CPRLOG)
      F = 10.0D0**( FCLOG/(1.0D0+X*X) )
      PCOR = PR*F/(1.0D0+PR)
      RF(70) = RFH70*PCOR                                                       
      RB(70) = EXP(+(-2.77264D+07)*TM3+( 3.43893D+05)*TM2                       
     &         +(-1.99235D+04)*TM1+( 3.00214D+01)+(-2.41747D-04)*TP1            
     &         +( 1.77849D-09)*TP2+( 2.94459D-12)*TP3)                          
      RB(70) = RB(70)*PCOR                                                      
      RFH71 = 6.0800D+12 * EXP(  0.270*ALOGT -  0.28003*RTR )                   
      RFL71 = 1.4000D+30 * EXP( -3.860*ALOGT -  3.32036*RTR )                   
      PR = RFL71*XM(71)/RFH71                                                   
      PRLOG = DLOG10(MAX(PR,SMALL))
      F4 =     7.82000D-01
      F5 =     2.07500D+02
      F6 =     2.66300D+03
      F7 =     6.09500D+03
      FC = (1.0D0-F4)*EXP(-T/F5)+F4*EXP(-T/F6)+EXP(-F7/T)
      FCLOG = DLOG10(MAX(FC,SMALL))
      CPRLOG = PRLOG - (0.4D0 + 0.67D0*FCLOG)
      X = CPRLOG/(0.75D0-1.27D0*FCLOG-0.14D0*CPRLOG)
      F = 10.0D0**( FCLOG/(1.0D0+X*X) )
      PCOR = PR*F/(1.0D0+PR)
      RF(71) = RFH71*PCOR                                                       
      RB(71) = EXP(+(-2.70894D+07)*TM3+( 3.35028D+05)*TM2                       
     &         +(-5.72335D+04)*TM1+( 3.75254D+01)+(-1.11866D-04)*TP1            
     &         +(-5.55496D-08)*TP2+( 8.81053D-12)*TP3)                          
      RB(71) = RB(71)*PCOR                                                      
      RFH73 = 1.0800D+12 * EXP(  0.454*ALOGT -  1.82020*RTR )                   
      RFL73 = 1.2000D+42 * EXP( -7.620*ALOGT -  6.97076*RTR )                   
      PR = RFL73*XM(73)/RFH73                                                   
      PRLOG = DLOG10(MAX(PR,SMALL))
      F4 =     9.75300D-01
      F5 =     2.10000D+02
      F6 =     9.84000D+02
      F7 =     4.37400D+03
      FC = (1.0D0-F4)*EXP(-T/F5)+F4*EXP(-T/F6)+EXP(-F7/T)
      FCLOG = DLOG10(MAX(FC,SMALL))
      CPRLOG = PRLOG - (0.4D0 + 0.67D0*FCLOG)
      X = CPRLOG/(0.75D0-1.27D0*FCLOG-0.14D0*CPRLOG)
      F = 10.0D0**( FCLOG/(1.0D0+X*X) )
      PCOR = PR*F/(1.0D0+PR)
      RF(73) = RFH73*PCOR                                                       
      RB(73) = EXP(+(-2.05168D+07)*TM3+( 2.25656D+05)*TM2                       
     &         +(-1.98366D+04)*TM1+( 3.10392D+01)+( 1.89973D-04)*TP1            
     &         +(-9.10833D-08)*TP2+( 1.11107D-11)*TP3)                          
      RB(73) = RB(73)*PCOR                                                      
      RFH75 = 5.2100D+17 * EXP( -0.990*ALOGT -  1.58017*RTR )                   
      RFL75 = 1.9900D+41 * EXP( -7.080*ALOGT -  6.68573*RTR )                   
      PR = RFL75*XM(75)/RFH75                                                   
      PRLOG = DLOG10(MAX(PR,SMALL))
      F4 =     8.42200D-01
      F5 =     1.25000D+02
      F6 =     2.21900D+03
      F7 =     6.88200D+03
      FC = (1.0D0-F4)*EXP(-T/F5)+F4*EXP(-T/F6)+EXP(-F7/T)
      FCLOG = DLOG10(MAX(FC,SMALL))
      CPRLOG = PRLOG - (0.4D0 + 0.67D0*FCLOG)
      X = CPRLOG/(0.75D0-1.27D0*FCLOG-0.14D0*CPRLOG)
      F = 10.0D0**( FCLOG/(1.0D0+X*X) )
      PCOR = PR*F/(1.0D0+PR)
      RF(75) = RFH75*PCOR                                                       
      RB(75) = EXP(+(-1.40102D+07)*TM3+( 1.66802D+05)*TM2                       
     &         +(-5.15163D+04)*TM1+( 4.04437D+01)+(-1.04942D-03)*TP1            
     &         +( 1.07581D-07)*TP2+(-5.45529D-12)*TP3)                          
      RB(75) = RB(75)*PCOR                                                      
      RFH82 = 4.3000D+07 * EXP(  1.500*ALOGT - 79.60864*RTR )                   
      RFL82 = 5.0700D+27 * EXP( -3.420*ALOGT - 84.35915*RTR )                   
      PR = RFL82*XM(82)/RFH82                                                   
      PRLOG = DLOG10(MAX(PR,SMALL))
      F4 =     9.32000D-01
      F5 =     1.97000D+02
      F6 =     1.54000D+03
      F7 =     1.03000D+04
      FC = (1.0D0-F4)*EXP(-T/F5)+F4*EXP(-T/F6)+EXP(-F7/T)
      FCLOG = DLOG10(MAX(FC,SMALL))
      CPRLOG = PRLOG - (0.4D0 + 0.67D0*FCLOG)
      X = CPRLOG/(0.75D0-1.27D0*FCLOG-0.14D0*CPRLOG)
      F = 10.0D0**( FCLOG/(1.0D0+X*X) )
      PCOR = PR*F/(1.0D0+PR)
      RF(82) = RFH82*PCOR                                                       
      RB(82) = EXP(+(-4.92553D+07)*TM3+( 6.04826D+05)*TM2                       
     &         +(-4.23911D+04)*TM1+( 3.28651D+01)+(-2.19077D-05)*TP1            
     &         +(-5.04236D-08)*TP2+( 9.08484D-12)*TP3)                          
      RB(82) = RB(82)*PCOR                                                      
      RFH84 = 7.4000D+13 * EXP( -0.370*ALOGT                )                   
      RFL84 = 2.3000D+18 * EXP( -0.900*ALOGT +  1.70018*RTR )                   
      PR = RFL84*XM(84)/RFH84                                                   
      PRLOG = DLOG10(MAX(PR,SMALL))
      F4 =     7.34600D-01
      F5 =     9.40000D+01
      F6 =     1.75600D+03
      F7 =     5.18200D+03
      FC = (1.0D0-F4)*EXP(-T/F5)+F4*EXP(-T/F6)+EXP(-F7/T)
      FCLOG = DLOG10(MAX(FC,SMALL))
      CPRLOG = PRLOG - (0.4D0 + 0.67D0*FCLOG)
      X = CPRLOG/(0.75D0-1.27D0*FCLOG-0.14D0*CPRLOG)
      F = 10.0D0**( FCLOG/(1.0D0+X*X) )
      PCOR = PR*F/(1.0D0+PR)
      RF(84) = RFH84*PCOR                                                       
      RB(84) = EXP(+(-7.24444D+06)*TM3+( 1.46842D+05)*TM2                       
     &         +(-2.62562D+04)*TM1+( 3.62729D+01)+(-1.52775D-03)*TP1            
     &         +( 2.65089D-07)*TP2+(-2.28363D-11)*TP3)                          
      RB(84) = RB(84)*PCOR                                                      
      RFH94 = 6.3000D+13                                                        
      RFL94 = 2.7000D+38 * EXP( -6.300*ALOGT -  3.10034*RTR )                   
      PR = RFL94*XM(94)/RFH94                                                   
      PRLOG = DLOG10(MAX(PR,SMALL))
      F4 =     2.10500D-01
      F5 =     8.35000D+01
      F6 =     5.39800D+03
      F7 =     8.37000D+03
      FC = (1.0D0-F4)*EXP(-T/F5)+F4*EXP(-T/F6)+EXP(-F7/T)
      FCLOG = DLOG10(MAX(FC,SMALL))
      CPRLOG = PRLOG - (0.4D0 + 0.67D0*FCLOG)
      X = CPRLOG/(0.75D0-1.27D0*FCLOG-0.14D0*CPRLOG)
      F = 10.0D0**( FCLOG/(1.0D0+X*X) )
      PCOR = PR*F/(1.0D0+PR)
      RF(94) = RFH94*PCOR                                                       
      RB(94) = EXP(+(-3.00978D+07)*TM3+( 4.11002D+05)*TM2                       
     &         +(-4.81230D+04)*TM1+( 4.03773D+01)+(-1.36918D-03)*TP1            
     &         +( 2.20682D-07)*TP2+(-1.74721D-11)*TP3)                          
      RB(94) = RB(94)*PCOR                                                      
      RFH130 = 5.0000D+13                                                       
      RFL130 = 2.6900D+28 * EXP( -3.740*ALOGT -  1.93621*RTR )                  
      PR = RFL130*XM(130)/RFH130                                                
      PRLOG = DLOG10(MAX(PR,SMALL))
      F4 =     5.75700D-01
      F5 =     2.37000D+02
      F6 =     1.65200D+03
      F7 =     5.06900D+03
      FC = (1.0D0-F4)*EXP(-T/F5)+F4*EXP(-T/F6)+EXP(-F7/T)
      FCLOG = DLOG10(MAX(FC,SMALL))
      CPRLOG = PRLOG - (0.4D0 + 0.67D0*FCLOG)
      X = CPRLOG/(0.75D0-1.27D0*FCLOG-0.14D0*CPRLOG)
      F = 10.0D0**( FCLOG/(1.0D0+X*X) )
      PCOR = PR*F/(1.0D0+PR)
      RF(130) = RFH130*PCOR                                                     
      RB(130) = EXP(+( 1.13935D+07)*TM3+(-5.02156D+04)*TM2                      
     &          +(-3.70628D+04)*TM1+( 3.75999D+01)+(-1.46106D-03)*TP1           
     &          +( 2.97031D-07)*TP2+(-2.75444D-11)*TP3)                         
      RB(130) = RB(130)*PCOR                                                    
      RFH139 = 8.1000D+11 * EXP(  0.500*ALOGT -  4.51049*RTR )                  
      RFL139 = 2.6900D+33 * EXP( -5.110*ALOGT -  7.09577*RTR )                  
      PR = RFL139*XM(139)/RFH139                                                
      PRLOG = DLOG10(MAX(PR,SMALL))
      F4 =     5.90700D-01
      F5 =     2.75000D+02
      F6 =     1.22600D+03
      F7 =     5.18500D+03
      FC = (1.0D0-F4)*EXP(-T/F5)+F4*EXP(-T/F6)+EXP(-F7/T)
      FCLOG = DLOG10(MAX(FC,SMALL))
      CPRLOG = PRLOG - (0.4D0 + 0.67D0*FCLOG)
      X = CPRLOG/(0.75D0-1.27D0*FCLOG-0.14D0*CPRLOG)
      F = 10.0D0**( FCLOG/(1.0D0+X*X) )
      PCOR = PR*F/(1.0D0+PR)
      RF(139) = RFH139*PCOR                                                     
      RB(139) = EXP(+(-2.40706D+06)*TM3+( 1.08055D+05)*TM2                      
     &          +(-4.24931D+04)*TM1+( 3.94623D+01)+(-1.30138D-03)*TP1           
     &          +( 2.14502D-07)*TP2+(-1.81351D-11)*TP3)                         
      RB(139) = RB(139)*PCOR                                                    
      RFH145 = 2.0000D+13                                                       
      RFL145 = 2.7000D+38 * EXP( -6.300*ALOGT -  3.10034*RTR )                  
      PR = RFL145*XM(145)/RFH145                                                
      PRLOG = DLOG10(MAX(PR,SMALL))
      F4 =     1.50700D-01
      F5 =     1.34000D+02
      F6 =     2.38300D+03
      F7 =     7.26500D+03
      FC = (1.0D0-F4)*EXP(-T/F5)+F4*EXP(-T/F6)+EXP(-F7/T)
      FCLOG = DLOG10(MAX(FC,SMALL))
      CPRLOG = PRLOG - (0.4D0 + 0.67D0*FCLOG)
      X = CPRLOG/(0.75D0-1.27D0*FCLOG-0.14D0*CPRLOG)
      F = 10.0D0**( FCLOG/(1.0D0+X*X) )
      PCOR = PR*F/(1.0D0+PR)
      RF(145) = RFH145*PCOR                                                     
      RB(145) = EXP(+(-2.59165D+07)*TM3+( 3.79417D+05)*TM2                      
     &          +(-4.83024D+04)*TM1+( 3.94532D+01)+(-1.80081D-03)*TP1           
     &          +( 3.02933D-07)*TP2+(-2.44290D-11)*TP3)                         
      RB(145) = RB(145)*PCOR                                                    
      RFH156 = 2.1200D+16 * EXP( -0.970*ALOGT -  0.62007*RTR )                  
      RFL156 = 1.7700D+50 * EXP( -9.670*ALOGT -  6.22067*RTR )                  
      PR = RFL156*XM(156)/RFH156                                                
      PRLOG = DLOG10(MAX(PR,SMALL))
      F4 =     5.32500D-01
      F5 =     1.51000D+02
      F6 =     1.03800D+03
      F7 =     4.97000D+03
      FC = (1.0D0-F4)*EXP(-T/F5)+F4*EXP(-T/F6)+EXP(-F7/T)
      FCLOG = DLOG10(MAX(FC,SMALL))
      CPRLOG = PRLOG - (0.4D0 + 0.67D0*FCLOG)
      X = CPRLOG/(0.75D0-1.27D0*FCLOG-0.14D0*CPRLOG)
      F = 10.0D0**( FCLOG/(1.0D0+X*X) )
      PCOR = PR*F/(1.0D0+PR)
      RF(156) = RFH156*PCOR                                                     
      RB(156) = EXP(+(-1.78486D+07)*TM3+( 2.86650D+05)*TM2                      
     &          +(-4.66662D+04)*TM1+( 4.23664D+01)+(-2.56171D-03)*TP1           
     &          +( 4.26863D-07)*TP2+(-3.55820D-11)*TP3)                         
      RB(156) = RB(156)*PCOR                                                    
      RFH172 = 8.0000D+12 * EXP(  0.440*ALOGT - 88.77963*RTR )                  
      RFL172 = 7.0000D+50 * EXP( -9.310*ALOGT - 99.87083*RTR )                  
      PR = RFL172*XM(172)/RFH172                                                
      PRLOG = DLOG10(MAX(PR,SMALL))
      F4 =     7.34500D-01
      F5 =     1.80000D+02
      F6 =     1.03500D+03
      F7 =     5.41700D+03
      FC = (1.0D0-F4)*EXP(-T/F5)+F4*EXP(-T/F6)+EXP(-F7/T)
      FCLOG = DLOG10(MAX(FC,SMALL))
      CPRLOG = PRLOG - (0.4D0 + 0.67D0*FCLOG)
      X = CPRLOG/(0.75D0-1.27D0*FCLOG-0.14D0*CPRLOG)
      F = 10.0D0**( FCLOG/(1.0D0+X*X) )
      PCOR = PR*F/(1.0D0+PR)
      RF(172) = RFH172*PCOR                                                     
      RB(172) = EXP(+( 4.29184D+07)*TM3+(-5.41368D+05)*TM2                      
     &          +(-2.16709D+04)*TM1+( 2.66179D+01)+( 1.36414D-03)*TP1           
     &          +(-1.83619D-07)*TP2+( 1.07724D-11)*TP3)                         
      RB(172) = RB(172)*PCOR                                                    
      RFH183 = 1.3000D+11 * EXP(              - 59.62647*RTR )                  
      RFL183 = 6.2000D+14 * EXP(              - 56.10609*RTR )                  
      PR = RFL183*XM(183)/RFH183                                                
      PCOR = PR/(1.0D0+PR)
      RF(183) = RFH183*PCOR                                                     
      RB(183) = EXP(+(-2.24252D+06)*TM3+(-3.63053D+04)*TM2                      
     &          +(-9.75088D+03)*TM1+( 1.96077D+01)+( 1.06890D-03)*TP1           
     &          +(-1.76739D-07)*TP2+( 1.55739D-11)*TP3)                         
      RB(183) = RB(183)*PCOR                                                    
      RFH239 = 3.1000D+12 * EXP(  0.150*ALOGT                )                  
      RFL239 = 1.3000D+25 * EXP( -3.160*ALOGT -  0.74008*RTR )                  
      PR = RFL239*XM(239)/RFH239                                                
      PRLOG = DLOG10(MAX(PR,SMALL))
      F4 =     6.67000D-01
      F5 =     2.35000D+02
      F6 =     2.11700D+03
      F7 =     4.53600D+03
      FC = (1.0D0-F4)*EXP(-T/F5)+F4*EXP(-T/F6)+EXP(-F7/T)
      FCLOG = DLOG10(MAX(FC,SMALL))
      CPRLOG = PRLOG - (0.4D0 + 0.67D0*FCLOG)
      X = CPRLOG/(0.75D0-1.27D0*FCLOG-0.14D0*CPRLOG)
      F = 10.0D0**( FCLOG/(1.0D0+X*X) )
      PCOR = PR*F/(1.0D0+PR)
      RF(239) = RFH239*PCOR                                                     
      RB(239) = EXP(+( 1.11118D+07)*TM3+(-3.71382D+04)*TM2                      
     &          +(-1.62466D+04)*TM1+( 3.48797D+01)+(-1.63042D-03)*TP1           
     &          +( 3.47564D-07)*TP2+(-3.28633D-11)*TP3)                         
      RB(239) = RB(239)*PCOR                                                    

      RETURN
      END

      SUBROUTINE STEADYE( ABV, DEN, RF, RB, XM, ADJ, CONMX, SMALL,
     &                    LBOUND, LITER, XH2, XH, XO, XO2, XOH, XH2O,           
     &                    XHO2, XH2O2, XC, XCH, XCH2, XCH2S, XCH3, XCH4,        
     &                    XCO, XCO2, XHCO, XCH2O, XCH2OH, XCH3O, XCH3OH,        
     &                    XC2H, XC2H2, XC2H3, XC2H4, XC2H5, XC2H6,              
     &                    XHCCO, XCH2CO, XHCCOH, XN, XNH, XNH2, XNH3,           
     &                    XNNH, XNO, XNO2, XN2O, XHNO, XCN, XHCN, XH2CN,        
     &                    XHCNN, XHCNO, XHOCN, XHNCO, XNCO, XN2 )               

C*****double precision
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END double precision
C*****single precision
C       IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END single precision

      DIMENSION ABV(*), DEN(*), RF(*), RB(*), XM(*)
      LOGICAL LBOUND(*), LITER

C     STEADY-STATE EXPRESSION FOR NNH       

      ABV(1) = RB(209)*XCH4*XN2 +RB(205)*XOH*XN2 +RB(204)*XHO2*XN2 +            
     &         RB(208)*XH2O*XN2 +RB(206)*XNH*XNO +RB(203)*XH*XN2*XM(203)        
     &         +RB(207)*XH2*XN2 +RB(202)*XH*XN2                                 
      DEN(1) = RF(209)*XCH3 +RF(205)*XO +RF(204)*XO2 +RF(208)*XOH +             
     &         RF(206)*XO +RF(203)*XM(203) +RF(207)*XH +RF(202)                 
      IF( LITER .OR. .NOT. LBOUND(1) ) THEN                                     
      IF(DEN(1) .LT. 1.0) DEN(1) = MAX( ADJ*ABV(1), DEN(1), SMALL )             
      VOLD = XNNH                                                               
      XNNH = ABV(1)/DEN(1)                                                      
      DIFF = ABS( (XNNH-VOLD)/MAX(XNNH,VOLD,SMALL) )                            
      CONMX = MAX( CONMX, DIFF )
      ENDIF

C     STEADY-STATE EXPRESSION FOR C2H       

      ABV(2) = RB(69)*XC2H2 +RF(22)*XO*XC2H2 +RB(170)*XH*XC2H2 +RB(105)         
     &         *XH*XHCCO +RB(20)*XCH*XCO +RF(122)*XC*XCH2 +RF(108)*XOH          
     &         *XC2H2 +RB(169)*XCO*XHCO                                         
      DEN(2) = RF(69)*XH +RB(22)*XOH +RF(170)*XH2 +RF(105)*XOH +RF(20)          
     &         *XO +RB(122)*XH +RB(108)*XH2O +RF(169)*XO2                       
      IF( LITER .OR. .NOT. LBOUND(2) ) THEN                                     
      IF(DEN(2) .LT. 1.0) DEN(2) = MAX( ADJ*ABV(2), DEN(2), SMALL )             
      VOLD = XC2H                                                               
      XC2H = ABV(2)/DEN(2)                                                      
      DIFF = ABS( (XC2H-VOLD)/MAX(XC2H,VOLD,SMALL) )                            
      CONMX = MAX( CONMX, DIFF )
      ENDIF

C     STEADY-STATE EXPRESSION FOR CN        

      ABV(3) = RF(228)*XHCN*XM(228) +RF(237)*XC*XN2 +RF(231)*XO*XHCN +          
     &         RB(219)*XH*XHCN +RF(242)*XC*XNO +RB(218)*XO*XNCO +RB(215)        
     &         *XCO*XN +RB(216)*XH*XNCO +RB(217)*XOH*XHCN                       
      DEN(3) = RB(228)*XH*XM(228) +RB(237)*XN +RB(231)*XOH +RF(219)*XH2         
     &         +RB(242)*XO +RF(218)*XO2 +RF(215)*XO +RF(216)*XOH +              
     &         RF(217)*XH2O                                                     
      IF( LITER .OR. .NOT. LBOUND(3) ) THEN                                     
      IF(DEN(3) .LT. 1.0) DEN(3) = MAX( ADJ*ABV(3), DEN(3), SMALL )             
      VOLD = XCN                                                                
      XCN = ABV(3)/DEN(3)                                                       
      DIFF = ABS( (XCN-VOLD)/MAX(XCN,VOLD,SMALL) )                              
      CONMX = MAX( CONMX, DIFF )
      ENDIF

C     STEADY-STATE EXPRESSION FOR HOCN      

      ABV(4) = RB(271)*XHNCO*XH +RF(232)*XOH*XHCN                               
      DEN(4) = RF(271)*XH +RB(232)*XH                                           
      IF( LITER .OR. .NOT. LBOUND(4) ) THEN                                     
      IF(DEN(4) .LT. 1.0) DEN(4) = MAX( ADJ*ABV(4), DEN(4), SMALL )             
      VOLD = XHOCN                                                              
      XHOCN = ABV(4)/DEN(4)                                                     
      DIFF = ABS( (XHOCN-VOLD)/MAX(XHOCN,VOLD,SMALL) )                          
      CONMX = MAX( CONMX, DIFF )
      ENDIF

C     STEADY-STATE EXPRESSION FOR HNO       

      ABV(5) = RF(261)*XO*XHNCO +RB(214)*XHO2*XNO +RF(192)*XO2*XNH +            
     &         RB(211)*XOH*XNO +RF(195)*XH2O*XNH +RF(199)*XO*XNH2 +             
     &         RF(190)*XOH*XNH +RB(213)*XH2O*XNO +RB(212)*XH2*XNO +             
     &         RF(210)*XH*XNO*XM(210)                                           
      DEN(5) = RB(261)*XCO +RF(214)*XO2 +RB(192)*XO +RF(211)*XO +RB(195)        
     &         *XH2 +RB(199)*XH +RB(190)*XH +RF(213)*XOH +RF(212)*XH +          
     &         RB(210)*XM(210)                                                  
      IF( LITER .OR. .NOT. LBOUND(5) ) THEN                                     
      IF(DEN(5) .LT. 1.0) DEN(5) = MAX( ADJ*ABV(5), DEN(5), SMALL )             
      VOLD = XHNO                                                               
      XHNO = ABV(5)/DEN(5)                                                      
      DIFF = ABS( (XHNO-VOLD)/MAX(XHNO,VOLD,SMALL) )                            
      CONMX = MAX( CONMX, DIFF )
      ENDIF

C     STEADY-STATE EXPRESSION FOR HCCOH     

      ABV(6) = RF(107)*XOH*XC2H2 +RB(81)*XCH2CO*XH                              
      DEN(6) = RB(107)*XH +RF(81)*XH                                            
      IF( LITER .OR. .NOT. LBOUND(6) ) THEN                                     
      IF(DEN(6) .LT. 1.0) DEN(6) = MAX( ADJ*ABV(6), DEN(6), SMALL )             
      VOLD = XHCCOH                                                             
      XHCCOH = ABV(6)/DEN(6)                                                    
      DIFF = ABS( (XHCCOH-VOLD)/MAX(XHCCOH,VOLD,SMALL) )                        
      CONMX = MAX( CONMX, DIFF )
      ENDIF

C     STEADY-STATE EXPRESSION FOR C2H5      

      ABV(7) = RF(163)*XCH3*XC2H6 +RF(152)*XCH2S*XC2H6 +RB(173)*XHO2            
     &         *XC2H4 +RB(76)*XH2*XC2H4 +RB(75)*XC2H6 +RF(27)*XO*XC2H6 +        
     &         RF(157)*XCH3*XCH3 +RB(26)*XCH3*XCH2O +RF(112)*XOH*XC2H6 +        
     &         RF(77)*XH*XC2H6 +RF(73)*XH*XC2H4                                 
      DEN(7) = RB(163)*XCH4 +RB(152)*XCH3 +RF(173)*XO2 +RF(76)*XH +             
     &         RF(75)*XH +RB(27)*XOH +RB(157)*XH +RF(26)*XO +RB(112)            
     &         *XH2O +RB(77)*XH2 +RB(73)                                        
      IF( LITER .OR. .NOT. LBOUND(7) ) THEN                                     
      IF(DEN(7) .LT. 1.0) DEN(7) = MAX( ADJ*ABV(7), DEN(7), SMALL )             
      VOLD = XC2H5                                                              
      XC2H5 = ABV(7)/DEN(7)                                                     
      DIFF = ABS( (XC2H5-VOLD)/MAX(XC2H5,VOLD,SMALL) )                          
      CONMX = MAX( CONMX, DIFF )
      ENDIF

C     STEADY-STATE EXPRESSION FOR CH2S      

      ABV(8) = RB(241)*XNH*XHCN +RB(152)*XCH3*XC2H5 +RB(251)*XOH*XHCN +         
     &         RB(252)*XH*XHCNO +RB(250)*XH*XHNCO +RB(147)*XH*XC2H4 +           
     &         RF(66)*XH*XCH3O +RB(148)*XCH3*XCH3 +RB(9)*XH*XHCO +RB(8)         
     &         *XH2*XCO +RB(145)*XCH3OH +RF(78)*XH*XHCCO +RF(61)*XH             
     &         *XCH2OH +RB(93)*XH*XCH2O +RB(50)*XH2*XCH +RB(149)*XCH2           
     &         *XCO +RB(150)*XCH2*XCO2 +RB(143)*XH2O*XCO +RB(142)*XH*XOH        
     &         *XCO +RB(151)*XCO*XCH2O +RB(144)*XH*XCH3 +RB(146)*XCH2           
     &         *XH2O +RB(141)*XCH2*XN2 +RF(96)*XOH*XCH3                         
      DEN(8) = RF(241)*XN2 +RF(152)*XC2H6 +RF(251)*XNO +RF(252)*XNO +           
     &         RF(250)*XNO +RF(147)*XCH3 +RB(66)*XH2O +RF(148)*XCH4 +           
     &         RF(9)*XO +RF(8)*XO +RF(145)*XH2O +RB(78)*XCO +RB(61)*XH2O        
     &         +RF(93)*XOH +RF(50)*XH +RF(149)*XCO +RF(150)*XCO2 +              
     &         RF(143)*XO2 +RF(142)*XO2 +RF(151)*XCO2 +RF(144)*XH2 +            
     &         RF(146)*XH2O +RF(141)*XN2 +RB(96)*XH2O                           
      IF( LITER .OR. .NOT. LBOUND(8) ) THEN                                     
      IF(DEN(8) .LT. 1.0) DEN(8) = MAX( ADJ*ABV(8), DEN(8), SMALL )             
      VOLD = XCH2S                                                              
      XCH2S = ABV(8)/DEN(8)                                                     
      DIFF = ABS( (XCH2S-VOLD)/MAX(XCH2S,VOLD,SMALL) )                          
      CONMX = MAX( CONMX, DIFF )
      ENDIF

C     STEADY-STATE EXPRESSION FOR C2H3      

      ABV(9) = RF(162)*XCH3*XC2H4 +RF(140)*XCH2*XHCCO +RB(71)*XC2H4 +           
     &         RB(110)*XH2O*XC2H2 +RB(24)*XH*XCH2CO +RF(128)*XCH*XCH3 +         
     &         RB(171)*XHCO*XCH2O +RB(72)*XH2*XC2H2 +RF(111)*XOH*XC2H4 +        
     &         RF(74)*XH*XC2H4 +RF(70)*XH*XC2H2                                 
      DEN(9) = RB(162)*XCH4 +RB(140)*XCO +RF(71)*XH +RF(110)*XOH +RF(24)        
     &         *XO +RB(128)*XH +RF(171)*XO2 +RF(72)*XH +RB(111)*XH2O +          
     &         RB(74)*XH2 +RB(70)                                               
      IF( LITER .OR. .NOT. LBOUND(9) ) THEN                                     
      IF(DEN(9) .LT. 1.0) DEN(9) = MAX( ADJ*ABV(9), DEN(9), SMALL )             
      VOLD = XC2H3                                                              
      XC2H3 = ABV(9)/DEN(9)                                                     
      DIFF = ABS( (XC2H3-VOLD)/MAX(XC2H3,VOLD,SMALL) )                          
      CONMX = MAX( CONMX, DIFF )
      ENDIF

C     STEADY-STATE EXPRESSION FOR CH2OH     

      ABV(10) = RF(160)*XCH3*XCH3OH +RB(58)*XCH3OH +RF(63)*XCH3O*XH +           
     &          RB(16)*XOH*XCH2O +RB(101)*XH2O*XCH2O +RF(103)*XOH*XCH3OH        
     &          +RB(61)*XH2O*XCH2S +RF(18)*XO*XCH3OH +RB(59)*XH2*XCH2O +        
     &          RB(167)*XHO2*XCH2O +RF(67)*XH*XCH3OH +RB(60)*XOH*XCH3 +         
     &          RF(55)*XH*XCH2O                                                 
      DEN(10) = RB(160)*XCH4 +RF(58)*XH +RB(63)*XH +RF(16)*XO +RF(101)          
     &          *XOH +RB(103)*XH2O +RF(61)*XH +RB(18)*XOH +RF(59)*XH +          
     &          RF(167)*XO2 +RB(67)*XH2 +RF(60)*XH +RB(55)                      
      IF( LITER .OR. .NOT. LBOUND(10) ) THEN                                    
      IF(DEN(10) .LT. 1.0) DEN(10) = MAX( ADJ*ABV(10), DEN(10), SMALL )         
      VOLD = XCH2OH                                                             
      XCH2OH = ABV(10)/DEN(10)                                                  
      DIFF = ABS( (XCH2OH-VOLD)/MAX(XCH2OH,VOLD,SMALL) )                        
      CONMX = MAX( CONMX, DIFF )
      ENDIF

C     STEADY-STATE EXPRESSION FOR CH2       

      ABV(11) = RF(236)*XN*XH2CN +RB(240)*XNH*XHCN +RB(140)*XCO*XC2H3 +         
     &          RB(248)*XOH*XHCN +RB(122)*XH*XC2H +RB(249)*XH*XHCNO +           
     &          RB(127)*XH*XC2H2 +RF(259)*XH*XHCNN +RB(116)*XOH*XCH2O +         
     &          RB(247)*XH*XHNCO +RF(30)*XO*XCH2CO +RB(138)*XCH3*XCH3 +         
     &          RB(136)*XH2*XC2H2 +RB(139)*XCH2CO +RB(49)*XCH3 +RF(23)          
     &          *XO*XC2H2 +RB(137)*XH*XC2H4 +RB(135)*XH*XCH3 +RF(149)           
     &          *XCH2S*XCO +RF(150)*XCH2S*XCO2 +RB(92)*XH2O*XCH +RB(91)         
     &          *XH*XCH2O +RF(95)*XOH*XCH3 +RB(7)*XH*XHCO +RB(134)*XOH          
     &          *XHCO +RF(146)*XCH2S*XH2O +RF(141)*XCH2S*XN2 +RF(125)           
     &          *XH2*XCH                                                        
      DEN(11) = RB(236)*XN2 +RF(240)*XN2 +RF(140)*XHCCO +RF(248)*XNO +          
     &          RF(122)*XC +RF(249)*XNO +RF(127)*XCH +RB(259)*XN2 +             
     &          RF(116)*XHO2 +RF(247)*XNO +RB(30)*XCO2 +RF(138)*XCH4 +          
     &          RF(136)*XCH2 +RF(139)*XCO +RF(49)*XH +RB(23)*XCO +              
     &          RF(137)*XCH3 +RF(135)*XH2 +RB(149)*XCO +RB(150)*XCO2 +          
     &          RF(92)*XOH +RF(91)*XOH +RB(95)*XH2O +RF(7)*XO +RF(134)          
     &          *XO2 +RB(146)*XH2O +RB(141)*XN2 +RB(125)*XH                     
      IF( LITER .OR. .NOT. LBOUND(11) ) THEN                                    
      IF(DEN(11) .LT. 1.0) DEN(11) = MAX( ADJ*ABV(11), DEN(11), SMALL )         
      VOLD = XCH2                                                               
      ACH2 = RF(136)                                                            
      BCH2 = DEN(11) - ACH2*XCH2                                                
      ZCH2 = ABV(11)/BCH2                                                       
      IF( ACH2 .GT. SMALL ) THEN                                                
        B2P4AC = BCH2*BCH2 + 4.0*ACH2*ABV(11)                                   
        XCH2 = 0.5*( SQRT(B2P4AC) - BCH2 )/ACH2                                 
        IF( XCH2 .LT. SMALL ) XCH2 = ZCH2                                       
      ELSE
        XCH2 = ZCH2                                                             
      ENDIF
      DIFF = ABS( (XCH2-VOLD)/MAX(XCH2,VOLD,SMALL) )                            
      CONMX = MAX( CONMX, DIFF )
      ENDIF

C     STEADY-STATE EXPRESSION FOR HNCO      

      ABV(12) = RB(261)*XCO*XHNO +RB(267)*XCO*XNH*XM(267) +RB(262)*XOH          
     &          *XNCO +RB(264)*XH2*XNCO +RB(265)*XH2O*XNCO +RF(233)*XOH         
     &          *XHCN +RB(260)*XCO2*XNH +RB(266)*XCO2*XNH2 +RF(250)             
     &          *XCH2S*XNO +RF(268)*XHCNO*XH +RF(271)*XHOCN*XH +RF(247)         
     &          *XCH2*XNO +RB(263)*XCO*XNH2                                     
      DEN(12) = RF(261)*XO +RF(267)*XM(267) +RF(262)*XO +RF(264)*XH +           
     &          RF(265)*XOH +RB(233)*XH +RF(260)*XO +RF(266)*XOH +              
     &          RB(250)*XH +RB(268)*XH +RB(271)*XH +RB(247)*XH +RF(263)         
     &          *XH                                                             
      IF( LITER .OR. .NOT. LBOUND(12) ) THEN                                    
      IF(DEN(12) .LT. 1.0) DEN(12) = MAX( ADJ*ABV(12), DEN(12), SMALL )         
      VOLD = XHNCO                                                              
      XHNCO = ABV(12)/DEN(12)                                                   
      DIFF = ABS( (XHNCO-VOLD)/MAX(XHNCO,VOLD,SMALL) )                          
      CONMX = MAX( CONMX, DIFF )
      ENDIF

C     STEADY-STATE EXPRESSION FOR NH2       

      ABV(13) = RF(276)*XOH*XNH3 +RF(277)*XO*XNH3 +RF(275)*XH*XNH3 +            
     &          RF(270)*XH*XHCNO +RF(234)*XOH*XHCN +RB(198)*XOH*XNH +           
     &          RF(266)*XOH*XHNCO +RB(201)*XH2O*XNH +RB(200)*XH2*XNH +          
     &          RB(199)*XH*XHNO +RF(263)*XH*XHNCO                               
      DEN(13) = RB(276)*XH2O +RB(277)*XOH +RB(275)*XH2 +RB(270)*XCO +           
     &          RB(234)*XCO +RF(198)*XO +RB(266)*XCO2 +RF(201)*XOH +            
     &          RF(200)*XH +RF(199)*XO +RB(263)*XCO                             
      IF( LITER .OR. .NOT. LBOUND(13) ) THEN                                    
      IF(DEN(13) .LT. 1.0) DEN(13) = MAX( ADJ*ABV(13), DEN(13), SMALL )         
      VOLD = XNH2                                                               
      XNH2 = ABV(13)/DEN(13)                                                    
      DIFF = ABS( (XNH2-VOLD)/MAX(XNH2,VOLD,SMALL) )                            
      CONMX = MAX( CONMX, DIFF )
      ENDIF

C     STEADY-STATE EXPRESSION FOR CH        

      ABV(14) = RB(133)*XCO*XC2H2 +RF(20)*XO*XC2H +RB(127)*XH*XC2H2 +           
     &          RB(245)*XH*XNCO +RB(239)*XHCNN +RB(246)*XHCO*XN +RB(244)        
     &          *XO*XHCN +RB(130)*XHCCO +RB(238)*XN*XHCN +RB(128)*XH            
     &          *XC2H3 +RB(132)*XH*XCH2CO +RB(129)*XH*XC2H4 +RF(50)*XH          
     &          *XCH2S +RB(131)*XCO*XHCO +RB(6)*XH*XCO +RB(90)*XH*XHCO +        
     &          RF(92)*XOH*XCH2 +RB(48)*XH2*XC +RB(124)*XO*XHCO +RB(126)        
     &          *XH*XCH2O +RB(125)*XH*XCH2                                      
      DEN(14) = RF(133)*XHCCO +RB(20)*XCO +RF(127)*XCH2 +RF(245)*XNO +          
     &          RF(239)*XN2 +RF(246)*XNO +RF(244)*XNO +RF(130)*XCO +            
     &          RF(238)*XN2 +RF(128)*XCH3 +RF(132)*XCH2O +RF(129)*XCH4 +        
     &          RB(50)*XH2 +RF(131)*XCO2 +RF(6)*XO +RF(90)*XOH +RB(92)          
     &          *XH2O +RF(48)*XH +RF(124)*XO2 +RF(126)*XH2O +RF(125)*XH2        
      IF( LITER .OR. .NOT. LBOUND(14) ) THEN                                    
      IF(DEN(14) .LT. 1.0) DEN(14) = MAX( ADJ*ABV(14), DEN(14), SMALL )         
      VOLD = XCH                                                                
      XCH = ABV(14)/DEN(14)                                                     
      DIFF = ABS( (XCH-VOLD)/MAX(XCH,VOLD,SMALL) )                              
      CONMX = MAX( CONMX, DIFF )
      ENDIF

C     STEADY-STATE EXPRESSION FOR C         

      ABV(15) = RB(237)*XN*XCN +RB(122)*XH*XC2H +RB(242)*XO*XCN +RB(243)        
     &          *XCO*XN +RB(123)*XH*XC2H2 +RB(89)*XH*XCO +RB(121)*XO*XCO        
     &          +RF(48)*XH*XCH                                                  
      DEN(15) = RF(237)*XN2 +RF(122)*XCH2 +RF(242)*XNO +RF(243)*XNO +           
     &          RF(123)*XCH3 +RF(89)*XOH +RF(121)*XO2 +RB(48)*XH2               
      IF( LITER .OR. .NOT. LBOUND(15) ) THEN                                    
      IF(DEN(15) .LT. 1.0) DEN(15) = MAX( ADJ*ABV(15), DEN(15), SMALL )         
      VOLD = XC                                                                 
      XC = ABV(15)/DEN(15)                                                      
      DIFF = ABS( (XC-VOLD)/MAX(XC,VOLD,SMALL) )                                
      CONMX = MAX( CONMX, DIFF )
      ENDIF

C     STEADY-STATE EXPRESSION FOR N         

      ABV(16) = RB(194)*XH*XN2 +RB(223)*XCO*XN2 +RB(236)*XCH2*XN2 +             
     &          RB(274)*XH2*XHCN +RF(237)*XC*XN2 +RF(225)*XNCO*XM(225) +        
     &          RB(177)*XO*XNO +RB(176)*XO*XN2 +RB(273)*XH*XH2CN +              
     &          RF(243)*XC*XNO +RF(215)*XO*XCN +RF(246)*XCH*XNO +RF(191)        
     &          *XOH*XNH +RF(189)*XH*XNH +RF(238)*XCH*XN2 +RB(178)*XH           
     &          *XNO                                                            
      DEN(16) = RF(194)*XNH +RF(223)*XNCO +RF(236)*XH2CN +RF(274)*XCH3 +        
     &          RB(237)*XCN +RB(225)*XCO*XM(225) +RF(177)*XO2 +RF(176)          
     &          *XNO +RF(273)*XCH3 +RB(243)*XCO +RB(215)*XCO +RB(246)           
     &          *XHCO +RB(191)*XH2O +RB(189)*XH2 +RB(238)*XHCN +RF(178)         
     &          *XOH                                                            
      IF( LITER .OR. .NOT. LBOUND(16) ) THEN                                    
      IF(DEN(16) .LT. 1.0) DEN(16) = MAX( ADJ*ABV(16), DEN(16), SMALL )         
      VOLD = XN                                                                 
      XN = ABV(16)/DEN(16)                                                      
      DIFF = ABS( (XN-VOLD)/MAX(XN,VOLD,SMALL) )                                
      CONMX = MAX( CONMX, DIFF )
      ENDIF

C     STEADY-STATE EXPRESSION FOR H2CN      

      ABV(17) = RB(236)*XCH2*XN2 +RF(254)*XCH3*XNO +RF(273)*XCH3*XN +           
     &          RF(235)*XH*XHCN*XM(235)                                         
      DEN(17) = RF(236)*XN +RB(254)*XOH +RB(273)*XH +RB(235)*XM(235)            
      IF( LITER .OR. .NOT. LBOUND(17) ) THEN                                    
      IF(DEN(17) .LT. 1.0) DEN(17) = MAX( ADJ*ABV(17), DEN(17), SMALL )         
      VOLD = XH2CN                                                              
      XH2CN = ABV(17)/DEN(17)                                                   
      DIFF = ABS( (XH2CN-VOLD)/MAX(XH2CN,VOLD,SMALL) )                          
      CONMX = MAX( CONMX, DIFF )
      ENDIF

C     STEADY-STATE EXPRESSION FOR HCNN      

      ABV(18) = RB(256)*XNO*XHCN +RB(255)*XH*XCO*XN2 +RB(258)*XH*XHCO           
     &          *XN2 +RB(257)*XO*XHCO*XN2 +RB(259)*XCH2*XN2 +RF(239)*XCH        
     &          *XN2                                                            
      DEN(18) = RF(256)*XO +RF(255)*XO +RF(258)*XOH +RF(257)*XO2 +              
     &          RF(259)*XH +RB(239)                                             
      IF( LITER .OR. .NOT. LBOUND(18) ) THEN                                    
      IF(DEN(18) .LT. 1.0) DEN(18) = MAX( ADJ*ABV(18), DEN(18), SMALL )         
      VOLD = XHCNN                                                              
      XHCNN = ABV(18)/DEN(18)                                                   
      DIFF = ABS( (XHCNN-VOLD)/MAX(XHCNN,VOLD,SMALL) )                          
      CONMX = MAX( CONMX, DIFF )
      ENDIF

C     STEADY-STATE EXPRESSION FOR HCNO      

      ABV(19) = RB(270)*XCO*XNH2 +RB(269)*XOH*XHCN +RF(272)*XHCCO*XNO +         
     &          RF(252)*XCH2S*XNO +RF(249)*XCH2*XNO +RB(268)*XHNCO*XH           
      DEN(19) = RF(270)*XH +RF(269)*XH +RB(272)*XCO +RB(252)*XH +RB(249)        
     &          *XH +RF(268)*XH                                                 
      IF( LITER .OR. .NOT. LBOUND(19) ) THEN                                    
      IF(DEN(19) .LT. 1.0) DEN(19) = MAX( ADJ*ABV(19), DEN(19), SMALL )         
      VOLD = XHCNO                                                              
      XHCNO = ABV(19)/DEN(19)                                                   
      DIFF = ABS( (XHCNO-VOLD)/MAX(XHCNO,VOLD,SMALL) )                          
      CONMX = MAX( CONMX, DIFF )
      ENDIF

C     STEADY-STATE EXPRESSION FOR HCO       

      ABV(20) = RF(258)*XOH*XHCNN +RF(120)*XHO2*XCH2O +RF(257)*XO2*XHCNN        
     &          +RF(246)*XCH*XNO +RB(53)*XCH2O +RF(169)*XO2*XC2H +RF(32)        
     &          *XO2*XCH2O +RF(171)*XO2*XC2H3 +RF(159)*XCH3*XCH2O +RF(9)        
     &          *XO*XCH2S +RF(25)*XO*XC2H4 +RB(158)*XCH4*XCO +RF(131)           
     &          *XCH*XCO2 +RF(90)*XOH*XCH +RB(13)*XOH*XCO +RB(14)*XH            
     &          *XCO2 +RF(124)*XO2*XCH +RF(134)*XO2*XCH2 +RF(7)*XO*XCH2         
     &          +RB(166)*XHO2*XCO +RB(99)*XH2O*XCO +RF(15)*XO*XCH2O +           
     &          RB(54)*XH2*XCO +RB(165)*XH*XCO*XM(165) +RF(100)*XOH             
     &          *XCH2O +RF(57)*XH*XCH2O +RB(164)*XH*XCO*XH2O                    
      DEN(20) = RB(258)*XH*XN2 +RB(120)*XH2O2 +RB(257)*XO*XN2 +RB(246)          
     &          *XN +RF(53)*XH +RB(169)*XCO +RB(32)*XHO2 +RB(171)*XCH2O         
     &          +RB(159)*XCH4 +RB(9)*XH +RB(25)*XCH3 +RF(158)*XCH3 +            
     &          RB(131)*XCO +RB(90)*XH +RF(13)*XO +RF(14)*XO +RB(124)*XO        
     &          +RB(134)*XOH +RB(7)*XH +RF(166)*XO2 +RF(99)*XOH +RB(15)         
     &          *XOH +RF(54)*XH +RF(165)*XM(165) +RB(100)*XH2O +RB(57)          
     &          *XH2 +RF(164)*XH2O                                              
      IF( LITER .OR. .NOT. LBOUND(20) ) THEN                                    
      IF(DEN(20) .LT. 1.0) DEN(20) = MAX( ADJ*ABV(20), DEN(20), SMALL )         
      VOLD = XHCO                                                               
      XHCO = ABV(20)/DEN(20)                                                    
      DIFF = ABS( (XHCO-VOLD)/MAX(XHCO,VOLD,SMALL) )                            
      CONMX = MAX( CONMX, DIFF )
      ENDIF

C     STEADY-STATE EXPRESSION FOR CH3OH     

      ABV(21) = RF(62)*XH*XCH3O +RB(161)*XCH4*XCH3O +RB(160)*XCH4*XCH2OH        
     &          +RF(58)*XH*XCH2OH +RB(19)*XOH*XCH3O +RF(145)*XH2O*XCH2S         
     &          +RB(103)*XH2O*XCH2OH +RB(68)*XH2*XCH3O +RB(18)*XOH              
     &          *XCH2OH +RB(104)*XH2O*XCH3O +RB(67)*XH2*XCH2OH +RF(94)          
     &          *XOH*XCH3                                                       
      DEN(21) = RB(62) +RF(161)*XCH3 +RF(160)*XCH3 +RB(58) +RF(19)*XO +         
     &          RB(145) +RF(103)*XOH +RF(68)*XH +RF(18)*XO +RF(104)*XOH         
     &          +RF(67)*XH +RB(94)                                              
      IF( LITER .OR. .NOT. LBOUND(21) ) THEN                                    
      IF(DEN(21) .LT. 1.0) DEN(21) = MAX( ADJ*ABV(21), DEN(21), SMALL )         
      VOLD = XCH3OH                                                             
      XCH3OH = ABV(21)/DEN(21)                                                  
      DIFF = ABS( (XCH3OH-VOLD)/MAX(XCH3OH,VOLD,SMALL) )                        
      CONMX = MAX( CONMX, DIFF )
      ENDIF

C     STEADY-STATE EXPRESSION FOR CH3O      

      ABV(22) = RB(62)*XCH3OH +RF(161)*XCH3*XCH3OH +RB(63)*XCH2OH*XH +          
     &          RB(17)*XOH*XCH2O +RB(102)*XH2O*XCH2O +RB(66)*XH2O*XCH2S         
     &          +RB(65)*XOH*XCH3 +RF(153)*XO2*XCH3 +RB(64)*XH2*XCH2O +          
     &          RB(168)*XHO2*XCH2O +RF(19)*XO*XCH3OH +RF(118)*XHO2*XCH3         
     &          +RF(68)*XH*XCH3OH +RF(104)*XOH*XCH3OH +RF(56)*XH*XCH2O          
      DEN(22) = RF(62)*XH +RB(161)*XCH4 +RF(63)*XH +RF(17)*XO +RF(102)          
     &          *XOH +RF(66)*XH +RF(65)*XH +RB(153)*XO +RF(64)*XH +             
     &          RF(168)*XO2 +RB(19)*XOH +RB(118)*XOH +RB(68)*XH2 +              
     &          RB(104)*XH2O +RB(56)                                            
      IF( LITER .OR. .NOT. LBOUND(22) ) THEN                                    
      IF(DEN(22) .LT. 1.0) DEN(22) = MAX( ADJ*ABV(22), DEN(22), SMALL )         
      VOLD = XCH3O                                                              
      XCH3O = ABV(22)/DEN(22)                                                   
      DIFF = ABS( (XCH3O-VOLD)/MAX(XCH3O,VOLD,SMALL) )                          
      CONMX = MAX( CONMX, DIFF )
      ENDIF

C     STEADY-STATE EXPRESSION FOR O         

      ABV(23) = RB(261)*XCO*XHNO +RB(277)*XOH*XNH2 +RB(256)*XNO*XHCN +          
     &          RB(180)*XNO*XNO +RB(179)*XO2*XN2 +RB(262)*XOH*XNCO +            
     &          RF(192)*XO2*XNH +RB(186)*XO2*XNO +RB(255)*XH*XCO*XN2 +          
     &          RB(260)*XCO2*XNH +RB(198)*XOH*XNH +RB(22)*XOH*XC2H +            
     &          RF(177)*XO2*XN +RB(20)*XCH*XCO +RF(257)*XO2*XHCNN +             
     &          RB(231)*XOH*XCN +RF(218)*XO2*XCN +RF(242)*XC*XNO +              
     &          RF(176)*XN*XNO +RB(205)*XOH*XN2 +RB(230)*XCO*XNH +              
     &          RB(220)*XCO*XNO +RB(215)*XCO*XN +RB(199)*XH*XHNO +RB(29)        
     &          *XOH*XHCCO +RB(211)*XOH*XNO +RB(30)*XCH2*XCO2 +RF(31)           
     &          *XO2*XCO +RB(17)*XOH*XCH2O +RB(206)*XNH*XNO +RF(183)            
     &          *XN2O +RB(188)*XH*XNO +RF(244)*XCH*XNO +RB(229)*XH*XNCO         
     &          +RB(185)*XNO2*XM(185) +RB(27)*XOH*XC2H5 +RF(153)*XO2            
     &          *XCH3 +RB(26)*XCH3*XCH2O +RB(24)*XH*XCH2CO +RB(19)*XOH          
     &          *XCH3O +RB(16)*XOH*XCH2O +RB(28)*XH*XCO*XCO +RB(1)*XO2          
     &          *XM(1) +RB(9)*XH*XHCO +RB(8)*XH2*XCO +RB(21)*XH*XHCCO +         
     &          RB(23)*XCH2*XCO +RB(5)*XOH*XHO2 +RB(25)*XCH3*XHCO +RB(2)        
     &          *XOH*XM(2) +RB(18)*XOH*XCH2OH +RF(43)*XH*XHO2 +RB(6)*XH         
     &          *XCO +RB(12)*XCO2*XM(12) +RB(4)*XO2*XOH +RB(13)*XOH*XCO         
     &          +RB(14)*XH*XCO2 +RF(121)*XO2*XC +RF(124)*XO2*XCH +RB(7)         
     &          *XH*XHCO +RB(15)*XOH*XHCO +RF(85)*XOH*XOH +RB(11)*XOH           
     &          *XCH3 +RB(3)*XH*XOH +RB(10)*XH*XCH2O +RF(37)*XH*XO2             
      DEN(23) = RF(261)*XHNCO +RF(277)*XNH3 +RF(256)*XHCNN +RF(180)*XN2O        
     &          +RF(179)*XN2O +RF(262)*XHNCO +RB(192)*XHNO +RF(186)*XNO2        
     &          +RF(255)*XHCNN +RF(260)*XHNCO +RF(198)*XNH2 +RF(22)             
     &          *XC2H2 +RB(177)*XNO +RF(20)*XC2H +RB(257)*XHCO*XN2 +            
     &          RF(231)*XHCN +RB(218)*XNCO +RB(242)*XCN +RB(176)*XN2 +          
     &          RF(205)*XNNH +RF(230)*XHCN +RF(220)*XNCO +RF(215)*XCN +         
     &          RF(199)*XNH2 +RF(29)*XCH2CO +RF(211)*XHNO +RF(30)*XCH2CO        
     &          +RB(31)*XCO2 +RF(17)*XCH3O +RF(206)*XNNH +RB(183)*XN2 +         
     &          RF(188)*XNH +RB(244)*XHCN +RF(229)*XHCN +RF(185)*XNO            
     &          *XM(185) +RF(27)*XC2H6 +RB(153)*XCH3O +RF(26)*XC2H5 +           
     &          RF(24)*XC2H3 +RF(19)*XCH3OH +RF(16)*XCH2OH +RF(28)*XHCCO        
     &          +RF(1)*XO*XM(1) +RF(9)*XCH2S +RF(8)*XCH2S +RF(21)*XC2H2         
     &          +RF(23)*XC2H2 +RF(5)*XH2O2 +RF(25)*XC2H4 +RF(2)*XH*XM(2)        
     &          +RF(18)*XCH3OH +RB(43)*XH2O +RF(6)*XCH +RF(12)*XCO              
     &          *XM(12) +RF(4)*XHO2 +RF(13)*XHCO +RF(14)*XHCO +RB(121)          
     &          *XCO +RB(124)*XHCO +RF(7)*XCH2 +RF(15)*XCH2O +RB(85)            
     &          *XH2O +RF(11)*XCH4 +RF(3)*XH2 +RF(10)*XCH3 +RB(37)*XOH          
      IF( LITER .OR. .NOT. LBOUND(23) ) THEN                                    
      IF(DEN(23) .LT. 1.0) DEN(23) = MAX( ADJ*ABV(23), DEN(23), SMALL )         
      VOLD = XO                                                                 
      AO = RF(1)*XM(1)                                                          
      BO = DEN(23) - AO*XO                                                      
      ZO = ABV(23)/BO                                                           
      IF( AO .GT. SMALL ) THEN                                                  
        B2P4AC = BO*BO + 4.0*AO*ABV(23)                                         
        XO = 0.5*( SQRT(B2P4AC) - BO )/AO                                       
        IF( XO .LT. SMALL ) XO = ZO                                             
      ELSE
        XO = ZO                                                                 
      ENDIF
      DIFF = ABS( (XO-VOLD)/MAX(XO,VOLD,SMALL) )                                
      CONMX = MAX( CONMX, DIFF )
      ENDIF

C     STEADY-STATE EXPRESSION FOR N2O       

      ABV(24) = RF(226)*XNO*XNCO +RB(182)*XHO2*XN2 +RB(180)*XNO*XNO +           
     &          RB(179)*XO2*XN2 +RF(197)*XNH*XNO +RB(181)*XOH*XN2 +             
     &          RB(183)*XO*XN2                                                  
      DEN(24) = RB(226)*XCO +RF(182)*XOH +RF(180)*XO +RF(179)*XO +              
     &          RB(197)*XH +RF(181)*XH +RF(183)                                 
      IF( LITER .OR. .NOT. LBOUND(24) ) THEN                                    
      IF(DEN(24) .LT. 1.0) DEN(24) = MAX( ADJ*ABV(24), DEN(24), SMALL )         
      VOLD = XN2O                                                               
      XN2O = ABV(24)/DEN(24)                                                    
      DIFF = ABS( (XN2O-VOLD)/MAX(XN2O,VOLD,SMALL) )                            
      CONMX = MAX( CONMX, DIFF )
      ENDIF

C     STEADY-STATE EXPRESSION FOR HCCO      

      ABV(25) = RB(175)*XCO*XCO*XC2H2 +RB(133)*XCO*XC2H2 +RB(140)*XCO           
     &          *XC2H3 +RB(272)*XCO*XHCNO +RF(105)*XOH*XC2H +RF(29)*XO          
     &          *XCH2CO +RB(174)*XOH*XCO*XCO +RF(130)*XCH*XCO +RF(113)          
     &          *XOH*XCH2CO +RF(79)*XH*XCH2CO +RB(28)*XH*XCO*XCO +RF(21)        
     &          *XO*XC2H2 +RB(78)*XCH2S*XCO                                     
      DEN(25) = RF(175)*XHCCO +RF(133)*XCH +RF(140)*XCH2 +RF(272)*XNO +         
     &          RB(105)*XH +RB(29)*XOH +RF(174)*XO2 +RB(130) +RB(113)           
     &          *XH2O +RB(79)*XH2 +RF(28)*XO +RB(21)*XH +RF(78)*XH              
      IF( LITER .OR. .NOT. LBOUND(25) ) THEN                                    
      IF(DEN(25) .LT. 1.0) DEN(25) = MAX( ADJ*ABV(25), DEN(25), SMALL )         
      VOLD = XHCCO                                                              
      AHCCO = RF(175)                                                           
      BHCCO = DEN(25) - AHCCO*XHCCO                                             
      ZHCCO = ABV(25)/BHCCO                                                     
      IF( AHCCO .GT. SMALL ) THEN                                               
        B2P4AC = BHCCO*BHCCO + 4.0*AHCCO*ABV(25)                                
        XHCCO = 0.5*( SQRT(B2P4AC) - BHCCO )/AHCCO                              
        IF( XHCCO .LT. SMALL ) XHCCO = ZHCCO                                    
      ELSE
        XHCCO = ZHCCO                                                           
      ENDIF
      DIFF = ABS( (XHCCO-VOLD)/MAX(XHCCO,VOLD,SMALL) )                          
      CONMX = MAX( CONMX, DIFF )
      ENDIF

C     STEADY-STATE EXPRESSION FOR NCO       

      ABV(26) = RB(223)*XCO*XN2 +RB(227)*XCO2*XN2 +RB(226)*XCO*XN2O +           
     &          RB(224)*XCO2*XNO +RF(262)*XO*XHNCO +RF(264)*XH*XHNCO +          
     &          RF(265)*XOH*XHNCO +RB(225)*XCO*XN*XM(225) +RB(222)*XH           
     &          *XCO*XNO +RF(218)*XO2*XCN +RB(220)*XCO*XNO +RF(245)*XCH         
     &          *XNO +RF(216)*XOH*XCN +RF(229)*XO*XHCN +RB(221)*XCO*XNH         
      DEN(26) = RF(223)*XN +RF(227)*XNO +RF(226)*XNO +RF(224)*XO2 +             
     &          RB(262)*XOH +RB(264)*XH2 +RB(265)*XH2O +RF(225)*XM(225)         
     &          +RF(222)*XOH +RB(218)*XO +RF(220)*XO +RB(245)*XH +              
     &          RB(216)*XH +RB(229)*XH +RF(221)*XH                              
      IF( LITER .OR. .NOT. LBOUND(26) ) THEN                                    
      IF(DEN(26) .LT. 1.0) DEN(26) = MAX( ADJ*ABV(26), DEN(26), SMALL )         
      VOLD = XNCO                                                               
      XNCO = ABV(26)/DEN(26)                                                    
      DIFF = ABS( (XNCO-VOLD)/MAX(XNCO,VOLD,SMALL) )                            
      CONMX = MAX( CONMX, DIFF )
      ENDIF

C     STEADY-STATE EXPRESSION FOR NH        

      ABV(27) = RF(241)*XCH2S*XN2 +RB(194)*XH*XN2 +RF(240)*XCH2*XN2 +           
     &          RF(267)*XHNCO*XM(267) +RB(196)*XOH*XN2 +RB(193)*XOH*XNO         
     &          +RB(192)*XO*XHNO +RF(260)*XO*XHNCO +RF(198)*XO*XNH2 +           
     &          RB(197)*XH*XN2O +RF(201)*XOH*XNH2 +RF(200)*XH*XNH2 +            
     &          RB(195)*XH2*XHNO +RF(230)*XO*XHCN +RB(190)*XH*XHNO +            
     &          RB(191)*XH2O*XN +RF(206)*XO*XNNH +RB(188)*XH*XNO +              
     &          RB(189)*XH2*XN +RF(221)*XH*XNCO                                 
      DEN(27) = RB(241)*XHCN +RF(194)*XN +RB(240)*XHCN +RB(267)*XCO             
     &          *XM(267) +RF(196)*XNO +RF(193)*XO2 +RF(192)*XO2 +RB(260)        
     &          *XCO2 +RB(198)*XOH +RF(197)*XNO +RB(201)*XH2O +RB(200)          
     &          *XH2 +RF(195)*XH2O +RB(230)*XCO +RF(190)*XOH +RF(191)           
     &          *XOH +RB(206)*XNO +RF(188)*XO +RF(189)*XH +RB(221)*XCO          
      IF( LITER .OR. .NOT. LBOUND(27) ) THEN                                    
      IF(DEN(27) .LT. 1.0) DEN(27) = MAX( ADJ*ABV(27), DEN(27), SMALL )         
      VOLD = XNH                                                                
      XNH = ABV(27)/DEN(27)                                                     
      DIFF = ABS( (XNH-VOLD)/MAX(XNH,VOLD,SMALL) )                              
      CONMX = MAX( CONMX, DIFF )
      ENDIF

C     STEADY-STATE EXPRESSION FOR NO2       

      ABV(28) = RB(186)*XO2*XNO +RF(184)*XHO2*XNO +RF(185)*XO*XNO               
     &          *XM(185) +RB(187)*XOH*XNO                                       
      DEN(28) = RF(186)*XO +RB(184)*XOH +RB(185)*XM(185) +RF(187)*XH            
      IF( LITER .OR. .NOT. LBOUND(28) ) THEN                                    
      IF(DEN(28) .LT. 1.0) DEN(28) = MAX( ADJ*ABV(28), DEN(28), SMALL )         
      VOLD = XNO2                                                               
      XNO2 = ABV(28)/DEN(28)                                                    
      DIFF = ABS( (XNO2-VOLD)/MAX(XNO2,VOLD,SMALL) )                            
      CONMX = MAX( CONMX, DIFF )
      ENDIF

C     STEADY-STATE EXPRESSION FOR CH2CO     

      ABV(29) = RF(81)*XHCCOH*XH +RF(106)*XOH*XC2H2 +RB(29)*XOH*XHCCO +         
     &          RB(30)*XCH2*XCO2 +RF(24)*XO*XC2H3 +RB(113)*XH2O*XHCCO +         
     &          RB(80)*XCH3*XCO +RF(139)*XCH2*XCO +RB(79)*XH2*XHCCO +           
     &          RF(132)*XCH*XCH2O                                               
      DEN(29) = RB(81)*XH +RB(106)*XH +RF(29)*XO +RF(30)*XO +RB(24)*XH +        
     &          RF(113)*XOH +RF(80)*XH +RB(139) +RF(79)*XH +RB(132)*XH          
      IF( LITER .OR. .NOT. LBOUND(29) ) THEN                                    
      IF(DEN(29) .LT. 1.0) DEN(29) = MAX( ADJ*ABV(29), DEN(29), SMALL )         
      VOLD = XCH2CO                                                             
      XCH2CO = ABV(29)/DEN(29)                                                  
      DIFF = ABS( (XCH2CO-VOLD)/MAX(XCH2CO,VOLD,SMALL) )                        
      CONMX = MAX( CONMX, DIFF )
      ENDIF

      RETURN
      END

      SUBROUTINE UPVALUE( IOPT, VNEW, BIG, LBOUND, CONMAX, XNNH, XC2H,          
     &                    XCN, XHOCN, XHNO, XHCCOH, XC2H5, XCH2S, XC2H3,        
     &                    XCH2OH, XCH2, XHNCO, XNH2, XCH, XC, XN, XH2CN,        
     &                    XHCNN, XHCNO, XHCO, XCH3OH, XCH3O, XO, XN2O,          
     &                    XHCCO, XNCO, XNH, XNO2, XCH2CO )                      

C*****double precision
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END double precision
C*****single precision
C       IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END single precision

      DIMENSION VNEW(*)
      LOGICAL LBOUND(*)

      IF( IOPT .EQ. 1 ) THEN
        XNNH = VNEW(1)                                                          
        XC2H = VNEW(2)                                                          
        XCN = VNEW(3)                                                           
        XHOCN = VNEW(4)                                                         
        XHNO = VNEW(5)                                                          
        XHCCOH = VNEW(6)                                                        
        XC2H5 = VNEW(7)                                                         
        XCH2S = VNEW(8)                                                         
        XC2H3 = VNEW(9)                                                         
        XCH2OH = VNEW(10)                                                       
        XCH2 = VNEW(11)                                                         
        XHNCO = VNEW(12)                                                        
        XNH2 = VNEW(13)                                                         
        XCH = VNEW(14)                                                          
        XC = VNEW(15)                                                           
        XN = VNEW(16)                                                           
        XH2CN = VNEW(17)                                                        
        XHCNN = VNEW(18)                                                        
        XHCNO = VNEW(19)                                                        
        XHCO = VNEW(20)                                                         
        XCH3OH = VNEW(21)                                                       
        XCH3O = VNEW(22)                                                        
        XO = VNEW(23)                                                           
        XN2O = VNEW(24)                                                         
        XHCCO = VNEW(25)                                                        
        XNCO = VNEW(26)                                                         
        XNH = VNEW(27)                                                          
        XNO2 = VNEW(28)                                                         
        XCH2CO = VNEW(29)                                                       
      ELSEIF( IOPT .EQ. 2 ) THEN
        CALL EXBOUND( XNNH, VNEW(1), LBOUND(1), CONMAX, BIG )                   
        CALL EXBOUND( XC2H, VNEW(2), LBOUND(2), CONMAX, BIG )                   
        CALL EXBOUND( XCN, VNEW(3), LBOUND(3), CONMAX, BIG )                    
        CALL EXBOUND( XHOCN, VNEW(4), LBOUND(4), CONMAX, BIG )                  
        CALL EXBOUND( XHNO, VNEW(5), LBOUND(5), CONMAX, BIG )                   
        CALL EXBOUND( XHCCOH, VNEW(6), LBOUND(6), CONMAX, BIG )                 
        CALL EXBOUND( XC2H5, VNEW(7), LBOUND(7), CONMAX, BIG )                  
        CALL EXBOUND( XCH2S, VNEW(8), LBOUND(8), CONMAX, BIG )                  
        CALL EXBOUND( XC2H3, VNEW(9), LBOUND(9), CONMAX, BIG )                  
        CALL EXBOUND( XCH2OH, VNEW(10), LBOUND(10), CONMAX, BIG )               
        CALL EXBOUND( XCH2, VNEW(11), LBOUND(11), CONMAX, BIG )                 
        CALL EXBOUND( XHNCO, VNEW(12), LBOUND(12), CONMAX, BIG )                
        CALL EXBOUND( XNH2, VNEW(13), LBOUND(13), CONMAX, BIG )                 
        CALL EXBOUND( XCH, VNEW(14), LBOUND(14), CONMAX, BIG )                  
        CALL EXBOUND( XC, VNEW(15), LBOUND(15), CONMAX, BIG )                   
        CALL EXBOUND( XN, VNEW(16), LBOUND(16), CONMAX, BIG )                   
        CALL EXBOUND( XH2CN, VNEW(17), LBOUND(17), CONMAX, BIG )                
        CALL EXBOUND( XHCNN, VNEW(18), LBOUND(18), CONMAX, BIG )                
        CALL EXBOUND( XHCNO, VNEW(19), LBOUND(19), CONMAX, BIG )                
        CALL EXBOUND( XHCO, VNEW(20), LBOUND(20), CONMAX, BIG )                 
        CALL EXBOUND( XCH3OH, VNEW(21), LBOUND(21), CONMAX, BIG )               
        CALL EXBOUND( XCH3O, VNEW(22), LBOUND(22), CONMAX, BIG )                
        CALL EXBOUND( XO, VNEW(23), LBOUND(23), CONMAX, BIG )                   
        CALL EXBOUND( XN2O, VNEW(24), LBOUND(24), CONMAX, BIG )                 
        CALL EXBOUND( XHCCO, VNEW(25), LBOUND(25), CONMAX, BIG )                
        CALL EXBOUND( XNCO, VNEW(26), LBOUND(26), CONMAX, BIG )                 
        CALL EXBOUND( XNH, VNEW(27), LBOUND(27), CONMAX, BIG )                  
        CALL EXBOUND( XNO2, VNEW(28), LBOUND(28), CONMAX, BIG )                 
        CALL EXBOUND( XCH2CO, VNEW(29), LBOUND(29), CONMAX, BIG )               
      ELSEIF( IOPT .EQ. 3 ) THEN
        VNEW(1) = XNNH                                                          
        VNEW(2) = XC2H                                                          
        VNEW(3) = XCN                                                           
        VNEW(4) = XHOCN                                                         
        VNEW(5) = XHNO                                                          
        VNEW(6) = XHCCOH                                                        
        VNEW(7) = XC2H5                                                         
        VNEW(8) = XCH2S                                                         
        VNEW(9) = XC2H3                                                         
        VNEW(10) = XCH2OH                                                       
        VNEW(11) = XCH2                                                         
        VNEW(12) = XHNCO                                                        
        VNEW(13) = XNH2                                                         
        VNEW(14) = XCH                                                          
        VNEW(15) = XC                                                           
        VNEW(16) = XN                                                           
        VNEW(17) = XH2CN                                                        
        VNEW(18) = XHCNN                                                        
        VNEW(19) = XHCNO                                                        
        VNEW(20) = XHCO                                                         
        VNEW(21) = XCH3OH                                                       
        VNEW(22) = XCH3O                                                        
        VNEW(23) = XO                                                           
        VNEW(24) = XN2O                                                         
        VNEW(25) = XHCCO                                                        
        VNEW(26) = XNCO                                                         
        VNEW(27) = XNH                                                          
        VNEW(28) = XNO2                                                         
        VNEW(29) = XCH2CO                                                       
      ENDIF

      RETURN
      END

      SUBROUTINE EXBOUND( VOLD, B, LBOUND, CONMAX, BIG )

C*****double precision
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END double precision
C*****single precision
C       IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END single precision

      PARAMETER (VBOUND = 100.0)

      LOGICAL LBOUND
      AA = LOG(VOLD) - B
      IF( ABS(AA) .LT. VBOUND ) THEN

        IF( EXP(AA) .LT. BIG ) THEN
          DIFF = ABS( (EXP(AA)-VOLD)/MAX(EXP(AA),VOLD) )
          CONMAX = MAX( CONMAX, DIFF )
          VOLD = EXP(AA)
          LBOUND = .TRUE.
        ELSE
          LBOUND = .FALSE.
        ENDIF

      ELSE
        LBOUND = .FALSE.
      ENDIF

      RETURN
      END

      SUBROUTINE NETRATE( W, RF, RB, XM, XH2, XH, XO, XO2, XOH, XH2O,           
     &                    XHO2, XH2O2, XC, XCH, XCH2, XCH2S, XCH3, XCH4,        
     &                    XCO, XCO2, XHCO, XCH2O, XCH2OH, XCH3O, XCH3OH,        
     &                    XC2H, XC2H2, XC2H3, XC2H4, XC2H5, XC2H6,              
     &                    XHCCO, XCH2CO, XHCCOH, XN, XNH, XNH2, XNH3,           
     &                    XNNH, XNO, XNO2, XN2O, XHNO, XCN, XHCN, XH2CN,        
     &                    XHCNN, XHCNO, XHOCN, XHNCO, XNCO, XN2 )               

C*****double precision
        IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END double precision
C*****single precision
C       IMPLICIT REAL (A-H, O-Z), INTEGER (I-N)
C*****END single precision

      DIMENSION W(*), RF(*), RB(*), XM(*)

C   NET PRODUCTION RATES FOR SKELETAL MECHANSIM

      W(1) = RF(1)*XO*XO*XM(1) - RB(1)*XO2*XM(1)                                
      W(2) = RF(2)*XH*XO*XM(2) - RB(2)*XOH*XM(2)                                
      W(3) = RF(3)*XH2*XO - RB(3)*XH*XOH                                        
      W(4) = RF(4)*XO*XHO2 - RB(4)*XO2*XOH                                      
      W(5) = RF(5)*XO*XH2O2 - RB(5)*XOH*XHO2                                    
      W(6) = RF(6)*XO*XCH - RB(6)*XH*XCO                                        
      W(7) = RF(7)*XO*XCH2 - RB(7)*XH*XHCO                                      
      W(8) = RF(8)*XO*XCH2S - RB(8)*XH2*XCO                                     
      W(9) = RF(9)*XO*XCH2S - RB(9)*XH*XHCO                                     
      W(10) = RF(10)*XO*XCH3 - RB(10)*XH*XCH2O                                  
      W(11) = RF(11)*XO*XCH4 - RB(11)*XOH*XCH3                                  
      W(12) = RF(12)*XO*XCO*XM(12) - RB(12)*XCO2*XM(12)                         
      W(13) = RF(13)*XO*XHCO - RB(13)*XOH*XCO                                   
      W(14) = RF(14)*XO*XHCO - RB(14)*XH*XCO2                                   
      W(15) = RF(15)*XO*XCH2O - RB(15)*XOH*XHCO                                 
      W(16) = RF(16)*XO*XCH2OH - RB(16)*XOH*XCH2O                               
      W(17) = RF(17)*XO*XCH3O - RB(17)*XOH*XCH2O                                
      W(18) = RF(18)*XO*XCH3OH - RB(18)*XOH*XCH2OH                              
      W(19) = RF(19)*XO*XCH3OH - RB(19)*XOH*XCH3O                               
      W(20) = RF(20)*XO*XC2H - RB(20)*XCH*XCO                                   
      W(21) = RF(21)*XO*XC2H2 - RB(21)*XH*XHCCO                                 
      W(22) = RF(22)*XO*XC2H2 - RB(22)*XOH*XC2H                                 
      W(23) = RF(23)*XO*XC2H2 - RB(23)*XCH2*XCO                                 
      W(24) = RF(24)*XO*XC2H3 - RB(24)*XH*XCH2CO                                
      W(25) = RF(25)*XO*XC2H4 - RB(25)*XCH3*XHCO                                
      W(26) = RF(26)*XO*XC2H5 - RB(26)*XCH3*XCH2O                               
      W(27) = RF(27)*XO*XC2H6 - RB(27)*XOH*XC2H5                                
      W(28) = RF(28)*XO*XHCCO - RB(28)*XH*XCO*XCO                               
      W(29) = RF(29)*XO*XCH2CO - RB(29)*XOH*XHCCO                               
      W(30) = RF(30)*XO*XCH2CO - RB(30)*XCH2*XCO2                               
      W(31) = RF(31)*XO2*XCO - RB(31)*XO*XCO2                                   
      W(32) = RF(32)*XO2*XCH2O - RB(32)*XHO2*XHCO                               
      W(33) = RF(33)*XH*XO2*XM(33) - RB(33)*XHO2*XM(33)                         
      W(34) = RF(34)*XH*XO2*XO2 - RB(34)*XHO2*XO2                               
      W(35) = RF(35)*XH*XO2*XH2O - RB(35)*XHO2*XH2O                             
      W(36) = RF(36)*XH*XO2*XN2 - RB(36)*XHO2*XN2                               
      W(37) = RF(37)*XH*XO2 - RB(37)*XO*XOH                                     
      W(38) = RF(38)*XH*XH*XM(38) - RB(38)*XH2*XM(38)                           
      W(39) = RF(39)*XH*XH*XH2 - RB(39)*XH2*XH2                                 
      W(40) = RF(40)*XH*XH*XH2O - RB(40)*XH2*XH2O                               
      W(41) = RF(41)*XH*XH*XCO2 - RB(41)*XH2*XCO2                               
      W(42) = RF(42)*XH*XOH*XM(42) - RB(42)*XH2O*XM(42)                         
      W(43) = RF(43)*XH*XHO2 - RB(43)*XO*XH2O                                   
      W(44) = RF(44)*XH*XHO2 - RB(44)*XH2*XO2                                   
      W(45) = RF(45)*XH*XHO2 - RB(45)*XOH*XOH                                   
      W(46) = RF(46)*XH*XH2O2 - RB(46)*XH2*XHO2                                 
      W(47) = RF(47)*XH*XH2O2 - RB(47)*XOH*XH2O                                 
      W(48) = RF(48)*XH*XCH - RB(48)*XH2*XC                                     
      W(49) = RF(49)*XH*XCH2 - RB(49)*XCH3                                      
      W(50) = RF(50)*XH*XCH2S - RB(50)*XH2*XCH                                  
      W(51) = RF(51)*XH*XCH3 - RB(51)*XCH4                                      
      W(52) = RF(52)*XH*XCH4 - RB(52)*XH2*XCH3                                  
      W(53) = RF(53)*XH*XHCO - RB(53)*XCH2O                                     
      W(54) = RF(54)*XH*XHCO - RB(54)*XH2*XCO                                   
      W(55) = RF(55)*XH*XCH2O - RB(55)*XCH2OH                                   
      W(56) = RF(56)*XH*XCH2O - RB(56)*XCH3O                                    
      W(57) = RF(57)*XH*XCH2O - RB(57)*XH2*XHCO                                 
      W(58) = RF(58)*XH*XCH2OH - RB(58)*XCH3OH                                  
      W(59) = RF(59)*XH*XCH2OH - RB(59)*XH2*XCH2O                               
      W(60) = RF(60)*XH*XCH2OH - RB(60)*XOH*XCH3                                
      W(61) = RF(61)*XH*XCH2OH - RB(61)*XH2O*XCH2S                              
      W(62) = RF(62)*XH*XCH3O - RB(62)*XCH3OH                                   
      W(63) = RF(63)*XCH3O*XH - RB(63)*XCH2OH*XH                                
      W(64) = RF(64)*XH*XCH3O - RB(64)*XH2*XCH2O                                
      W(65) = RF(65)*XH*XCH3O - RB(65)*XOH*XCH3                                 
      W(66) = RF(66)*XH*XCH3O - RB(66)*XH2O*XCH2S                               
      W(67) = RF(67)*XH*XCH3OH - RB(67)*XH2*XCH2OH                              
      W(68) = RF(68)*XH*XCH3OH - RB(68)*XH2*XCH3O                               
      W(69) = RF(69)*XH*XC2H - RB(69)*XC2H2                                     
      W(70) = RF(70)*XH*XC2H2 - RB(70)*XC2H3                                    
      W(71) = RF(71)*XH*XC2H3 - RB(71)*XC2H4                                    
      W(72) = RF(72)*XH*XC2H3 - RB(72)*XH2*XC2H2                                
      W(73) = RF(73)*XH*XC2H4 - RB(73)*XC2H5                                    
      W(74) = RF(74)*XH*XC2H4 - RB(74)*XH2*XC2H3                                
      W(75) = RF(75)*XH*XC2H5 - RB(75)*XC2H6                                    
      W(76) = RF(76)*XH*XC2H5 - RB(76)*XH2*XC2H4                                
      W(77) = RF(77)*XH*XC2H6 - RB(77)*XH2*XC2H5                                
      W(78) = RF(78)*XH*XHCCO - RB(78)*XCH2S*XCO                                
      W(79) = RF(79)*XH*XCH2CO - RB(79)*XH2*XHCCO                               
      W(80) = RF(80)*XH*XCH2CO - RB(80)*XCH3*XCO                                
      W(81) = RF(81)*XHCCOH*XH - RB(81)*XCH2CO*XH                               
      W(82) = RF(82)*XH2*XCO - RB(82)*XCH2O                                     
      W(83) = RF(83)*XH2*XOH - RB(83)*XH*XH2O                                   
      W(84) = RF(84)*XOH*XOH - RB(84)*XH2O2                                     
      W(85) = RF(85)*XOH*XOH - RB(85)*XO*XH2O                                   
      W(86) = RF(86)*XOH*XHO2 - RB(86)*XO2*XH2O                                 
      W(87) = RF(87)*XOH*XH2O2 - RB(87)*XH2O*XHO2                               
      W(88) = RF(88)*XOH*XH2O2 - RB(88)*XH2O*XHO2                               
      W(89) = RF(89)*XOH*XC - RB(89)*XH*XCO                                     
      W(90) = RF(90)*XOH*XCH - RB(90)*XH*XHCO                                   
      W(91) = RF(91)*XOH*XCH2 - RB(91)*XH*XCH2O                                 
      W(92) = RF(92)*XOH*XCH2 - RB(92)*XH2O*XCH                                 
      W(93) = RF(93)*XOH*XCH2S - RB(93)*XH*XCH2O                                
      W(94) = RF(94)*XOH*XCH3 - RB(94)*XCH3OH                                   
      W(95) = RF(95)*XOH*XCH3 - RB(95)*XH2O*XCH2                                
      W(96) = RF(96)*XOH*XCH3 - RB(96)*XH2O*XCH2S                               
      W(97) = RF(97)*XOH*XCH4 - RB(97)*XH2O*XCH3                                
      W(98) = RF(98)*XOH*XCO - RB(98)*XH*XCO2                                   
      W(99) = RF(99)*XOH*XHCO - RB(99)*XH2O*XCO                                 
      W(100) = RF(100)*XOH*XCH2O - RB(100)*XH2O*XHCO                            
      W(101) = RF(101)*XOH*XCH2OH - RB(101)*XH2O*XCH2O                          
      W(102) = RF(102)*XOH*XCH3O - RB(102)*XH2O*XCH2O                           
      W(103) = RF(103)*XOH*XCH3OH - RB(103)*XH2O*XCH2OH                         
      W(104) = RF(104)*XOH*XCH3OH - RB(104)*XH2O*XCH3O                          
      W(105) = RF(105)*XOH*XC2H - RB(105)*XH*XHCCO                              
      W(106) = RF(106)*XOH*XC2H2 - RB(106)*XH*XCH2CO                            
      W(107) = RF(107)*XOH*XC2H2 - RB(107)*XH*XHCCOH                            
      W(108) = RF(108)*XOH*XC2H2 - RB(108)*XH2O*XC2H                            
      W(109) = RF(109)*XOH*XC2H2 - RB(109)*XCH3*XCO                             
      W(110) = RF(110)*XOH*XC2H3 - RB(110)*XH2O*XC2H2                           
      W(111) = RF(111)*XOH*XC2H4 - RB(111)*XH2O*XC2H3                           
      W(112) = RF(112)*XOH*XC2H6 - RB(112)*XH2O*XC2H5                           
      W(113) = RF(113)*XOH*XCH2CO - RB(113)*XH2O*XHCCO                          
      W(114) = RF(114)*XHO2*XHO2 - RB(114)*XO2*XH2O2                            
      W(115) = RF(115)*XHO2*XHO2 - RB(115)*XO2*XH2O2                            
      W(116) = RF(116)*XHO2*XCH2 - RB(116)*XOH*XCH2O                            
      W(117) = RF(117)*XHO2*XCH3 - RB(117)*XO2*XCH4                             
      W(118) = RF(118)*XHO2*XCH3 - RB(118)*XOH*XCH3O                            
      W(119) = RF(119)*XHO2*XCO - RB(119)*XOH*XCO2                              
      W(120) = RF(120)*XHO2*XCH2O - RB(120)*XH2O2*XHCO                          
      W(121) = RF(121)*XO2*XC - RB(121)*XO*XCO                                  
      W(122) = RF(122)*XC*XCH2 - RB(122)*XH*XC2H                                
      W(123) = RF(123)*XC*XCH3 - RB(123)*XH*XC2H2                               
      W(124) = RF(124)*XO2*XCH - RB(124)*XO*XHCO                                
      W(125) = RF(125)*XH2*XCH - RB(125)*XH*XCH2                                
      W(126) = RF(126)*XH2O*XCH - RB(126)*XH*XCH2O                              
      W(127) = RF(127)*XCH*XCH2 - RB(127)*XH*XC2H2                              
      W(128) = RF(128)*XCH*XCH3 - RB(128)*XH*XC2H3                              
      W(129) = RF(129)*XCH*XCH4 - RB(129)*XH*XC2H4                              
      W(130) = RF(130)*XCH*XCO - RB(130)*XHCCO                                  
      W(131) = RF(131)*XCH*XCO2 - RB(131)*XCO*XHCO                              
      W(132) = RF(132)*XCH*XCH2O - RB(132)*XH*XCH2CO                            
      W(133) = RF(133)*XCH*XHCCO - RB(133)*XCO*XC2H2                            
      W(134) = RF(134)*XO2*XCH2 - RB(134)*XOH*XHCO                              
      W(135) = RF(135)*XH2*XCH2 - RB(135)*XH*XCH3                               
      W(136) = RF(136)*XCH2*XCH2 - RB(136)*XH2*XC2H2                            
      W(137) = RF(137)*XCH2*XCH3 - RB(137)*XH*XC2H4                             
      W(138) = RF(138)*XCH2*XCH4 - RB(138)*XCH3*XCH3                            
      W(139) = RF(139)*XCH2*XCO - RB(139)*XCH2CO                                
      W(140) = RF(140)*XCH2*XHCCO - RB(140)*XCO*XC2H3                           
      W(141) = RF(141)*XCH2S*XN2 - RB(141)*XCH2*XN2                             
      W(142) = RF(142)*XO2*XCH2S - RB(142)*XH*XOH*XCO                           
      W(143) = RF(143)*XO2*XCH2S - RB(143)*XH2O*XCO                             
      W(144) = RF(144)*XH2*XCH2S - RB(144)*XH*XCH3                              
      W(145) = RF(145)*XH2O*XCH2S - RB(145)*XCH3OH                              
      W(146) = RF(146)*XCH2S*XH2O - RB(146)*XCH2*XH2O                           
      W(147) = RF(147)*XCH2S*XCH3 - RB(147)*XH*XC2H4                            
      W(148) = RF(148)*XCH2S*XCH4 - RB(148)*XCH3*XCH3                           
      W(149) = RF(149)*XCH2S*XCO - RB(149)*XCH2*XCO                             
      W(150) = RF(150)*XCH2S*XCO2 - RB(150)*XCH2*XCO2                           
      W(151) = RF(151)*XCH2S*XCO2 - RB(151)*XCO*XCH2O                           
      W(152) = RF(152)*XCH2S*XC2H6 - RB(152)*XCH3*XC2H5                         
      W(153) = RF(153)*XO2*XCH3 - RB(153)*XO*XCH3O                              
      W(154) = RF(154)*XO2*XCH3 - RB(154)*XOH*XCH2O                             
      W(155) = RF(155)*XH2O2*XCH3 - RB(155)*XHO2*XCH4                           
      W(156) = RF(156)*XCH3*XCH3 - RB(156)*XC2H6                                
      W(157) = RF(157)*XCH3*XCH3 - RB(157)*XH*XC2H5                             
      W(158) = RF(158)*XCH3*XHCO - RB(158)*XCH4*XCO                             
      W(159) = RF(159)*XCH3*XCH2O - RB(159)*XCH4*XHCO                           
      W(160) = RF(160)*XCH3*XCH3OH - RB(160)*XCH4*XCH2OH                        
      W(161) = RF(161)*XCH3*XCH3OH - RB(161)*XCH4*XCH3O                         
      W(162) = RF(162)*XCH3*XC2H4 - RB(162)*XCH4*XC2H3                          
      W(163) = RF(163)*XCH3*XC2H6 - RB(163)*XCH4*XC2H5                          
      W(164) = RF(164)*XHCO*XH2O - RB(164)*XH*XCO*XH2O                          
      W(165) = RF(165)*XHCO*XM(165) - RB(165)*XH*XCO*XM(165)                    
      W(166) = RF(166)*XO2*XHCO - RB(166)*XHO2*XCO                              
      W(167) = RF(167)*XO2*XCH2OH - RB(167)*XHO2*XCH2O                          
      W(168) = RF(168)*XO2*XCH3O - RB(168)*XHO2*XCH2O                           
      W(169) = RF(169)*XO2*XC2H - RB(169)*XCO*XHCO                              
      W(170) = RF(170)*XH2*XC2H - RB(170)*XH*XC2H2                              
      W(171) = RF(171)*XO2*XC2H3 - RB(171)*XHCO*XCH2O                           
      W(172) = RF(172)*XC2H4 - RB(172)*XH2*XC2H2                                
      W(173) = RF(173)*XO2*XC2H5 - RB(173)*XHO2*XC2H4                           
      W(174) = RF(174)*XO2*XHCCO - RB(174)*XOH*XCO*XCO                          
      W(175) = RF(175)*XHCCO*XHCCO - RB(175)*XCO*XCO*XC2H2                      
      W(176) = RF(176)*XN*XNO - RB(176)*XO*XN2                                  
      W(177) = RF(177)*XO2*XN - RB(177)*XO*XNO                                  
      W(178) = RF(178)*XOH*XN - RB(178)*XH*XNO                                  
      W(179) = RF(179)*XO*XN2O - RB(179)*XO2*XN2                                
      W(180) = RF(180)*XO*XN2O - RB(180)*XNO*XNO                                
      W(181) = RF(181)*XH*XN2O - RB(181)*XOH*XN2                                
      W(182) = RF(182)*XOH*XN2O - RB(182)*XHO2*XN2                              
      W(183) = RF(183)*XN2O - RB(183)*XO*XN2                                    
      W(184) = RF(184)*XHO2*XNO - RB(184)*XOH*XNO2                              
      W(185) = RF(185)*XO*XNO*XM(185) - RB(185)*XNO2*XM(185)                    
      W(186) = RF(186)*XO*XNO2 - RB(186)*XO2*XNO                                
      W(187) = RF(187)*XH*XNO2 - RB(187)*XOH*XNO                                
      W(188) = RF(188)*XO*XNH - RB(188)*XH*XNO                                  
      W(189) = RF(189)*XH*XNH - RB(189)*XH2*XN                                  
      W(190) = RF(190)*XOH*XNH - RB(190)*XH*XHNO                                
      W(191) = RF(191)*XOH*XNH - RB(191)*XH2O*XN                                
      W(192) = RF(192)*XO2*XNH - RB(192)*XO*XHNO                                
      W(193) = RF(193)*XO2*XNH - RB(193)*XOH*XNO                                
      W(194) = RF(194)*XN*XNH - RB(194)*XH*XN2                                  
      W(195) = RF(195)*XH2O*XNH - RB(195)*XH2*XHNO                              
      W(196) = RF(196)*XNH*XNO - RB(196)*XOH*XN2                                
      W(197) = RF(197)*XNH*XNO - RB(197)*XH*XN2O                                
      W(198) = RF(198)*XO*XNH2 - RB(198)*XOH*XNH                                
      W(199) = RF(199)*XO*XNH2 - RB(199)*XH*XHNO                                
      W(200) = RF(200)*XH*XNH2 - RB(200)*XH2*XNH                                
      W(201) = RF(201)*XOH*XNH2 - RB(201)*XH2O*XNH                              
      W(202) = RF(202)*XNNH - RB(202)*XH*XN2                                    
      W(203) = RF(203)*XNNH*XM(203) - RB(203)*XH*XN2*XM(203)                    
      W(204) = RF(204)*XO2*XNNH - RB(204)*XHO2*XN2                              
      W(205) = RF(205)*XO*XNNH - RB(205)*XOH*XN2                                
      W(206) = RF(206)*XO*XNNH - RB(206)*XNH*XNO                                
      W(207) = RF(207)*XH*XNNH - RB(207)*XH2*XN2                                
      W(208) = RF(208)*XOH*XNNH - RB(208)*XH2O*XN2                              
      W(209) = RF(209)*XCH3*XNNH - RB(209)*XCH4*XN2                             
      W(210) = RF(210)*XH*XNO*XM(210) - RB(210)*XHNO*XM(210)                    
      W(211) = RF(211)*XO*XHNO - RB(211)*XOH*XNO                                
      W(212) = RF(212)*XH*XHNO - RB(212)*XH2*XNO                                
      W(213) = RF(213)*XOH*XHNO - RB(213)*XH2O*XNO                              
      W(214) = RF(214)*XO2*XHNO - RB(214)*XHO2*XNO                              
      W(215) = RF(215)*XO*XCN - RB(215)*XCO*XN                                  
      W(216) = RF(216)*XOH*XCN - RB(216)*XH*XNCO                                
      W(217) = RF(217)*XH2O*XCN - RB(217)*XOH*XHCN                              
      W(218) = RF(218)*XO2*XCN - RB(218)*XO*XNCO                                
      W(219) = RF(219)*XH2*XCN - RB(219)*XH*XHCN                                
      W(220) = RF(220)*XO*XNCO - RB(220)*XCO*XNO                                
      W(221) = RF(221)*XH*XNCO - RB(221)*XCO*XNH                                
      W(222) = RF(222)*XOH*XNCO - RB(222)*XH*XCO*XNO                            
      W(223) = RF(223)*XN*XNCO - RB(223)*XCO*XN2                                
      W(224) = RF(224)*XO2*XNCO - RB(224)*XCO2*XNO                              
      W(225) = RF(225)*XNCO*XM(225) - RB(225)*XCO*XN*XM(225)                    
      W(226) = RF(226)*XNO*XNCO - RB(226)*XCO*XN2O                              
      W(227) = RF(227)*XNO*XNCO - RB(227)*XCO2*XN2                              
      W(228) = RF(228)*XHCN*XM(228) - RB(228)*XH*XCN*XM(228)                    
      W(229) = RF(229)*XO*XHCN - RB(229)*XH*XNCO                                
      W(230) = RF(230)*XO*XHCN - RB(230)*XCO*XNH                                
      W(231) = RF(231)*XO*XHCN - RB(231)*XOH*XCN                                
      W(232) = RF(232)*XOH*XHCN - RB(232)*XH*XHOCN                              
      W(233) = RF(233)*XOH*XHCN - RB(233)*XH*XHNCO                              
      W(234) = RF(234)*XOH*XHCN - RB(234)*XCO*XNH2                              
      W(235) = RF(235)*XH*XHCN*XM(235) - RB(235)*XH2CN*XM(235)                  
      W(236) = RF(236)*XN*XH2CN - RB(236)*XCH2*XN2                              
      W(237) = RF(237)*XC*XN2 - RB(237)*XN*XCN                                  
      W(238) = RF(238)*XCH*XN2 - RB(238)*XN*XHCN                                
      W(239) = RF(239)*XCH*XN2 - RB(239)*XHCNN                                  
      W(240) = RF(240)*XCH2*XN2 - RB(240)*XNH*XHCN                              
      W(241) = RF(241)*XCH2S*XN2 - RB(241)*XNH*XHCN                             
      W(242) = RF(242)*XC*XNO - RB(242)*XO*XCN                                  
      W(243) = RF(243)*XC*XNO - RB(243)*XCO*XN                                  
      W(244) = RF(244)*XCH*XNO - RB(244)*XO*XHCN                                
      W(245) = RF(245)*XCH*XNO - RB(245)*XH*XNCO                                
      W(246) = RF(246)*XCH*XNO - RB(246)*XHCO*XN                                
      W(247) = RF(247)*XCH2*XNO - RB(247)*XH*XHNCO                              
      W(248) = RF(248)*XCH2*XNO - RB(248)*XOH*XHCN                              
      W(249) = RF(249)*XCH2*XNO - RB(249)*XH*XHCNO                              
      W(250) = RF(250)*XCH2S*XNO - RB(250)*XH*XHNCO                             
      W(251) = RF(251)*XCH2S*XNO - RB(251)*XOH*XHCN                             
      W(252) = RF(252)*XCH2S*XNO - RB(252)*XH*XHCNO                             
      W(253) = RF(253)*XCH3*XNO - RB(253)*XH2O*XHCN                             
      W(254) = RF(254)*XCH3*XNO - RB(254)*XOH*XH2CN                             
      W(255) = RF(255)*XO*XHCNN - RB(255)*XH*XCO*XN2                            
      W(256) = RF(256)*XO*XHCNN - RB(256)*XNO*XHCN                              
      W(257) = RF(257)*XO2*XHCNN - RB(257)*XO*XHCO*XN2                          
      W(258) = RF(258)*XOH*XHCNN - RB(258)*XH*XHCO*XN2                          
      W(259) = RF(259)*XH*XHCNN - RB(259)*XCH2*XN2                              
      W(260) = RF(260)*XO*XHNCO - RB(260)*XCO2*XNH                              
      W(261) = RF(261)*XO*XHNCO - RB(261)*XCO*XHNO                              
      W(262) = RF(262)*XO*XHNCO - RB(262)*XOH*XNCO                              
      W(263) = RF(263)*XH*XHNCO - RB(263)*XCO*XNH2                              
      W(264) = RF(264)*XH*XHNCO - RB(264)*XH2*XNCO                              
      W(265) = RF(265)*XOH*XHNCO - RB(265)*XH2O*XNCO                            
      W(266) = RF(266)*XOH*XHNCO - RB(266)*XCO2*XNH2                            
      W(267) = RF(267)*XHNCO*XM(267) - RB(267)*XCO*XNH*XM(267)                  
      W(268) = RF(268)*XHCNO*XH - RB(268)*XHNCO*XH                              
      W(269) = RF(269)*XH*XHCNO - RB(269)*XOH*XHCN                              
      W(270) = RF(270)*XH*XHCNO - RB(270)*XCO*XNH2                              
      W(271) = RF(271)*XHOCN*XH - RB(271)*XHNCO*XH                              
      W(272) = RF(272)*XHCCO*XNO - RB(272)*XCO*XHCNO                            
      W(273) = RF(273)*XCH3*XN - RB(273)*XH*XH2CN                               
      W(274) = RF(274)*XCH3*XN - RB(274)*XH2*XHCN                               
      W(275) = RF(275)*XH*XNH3 - RB(275)*XH2*XNH2                               
      W(276) = RF(276)*XOH*XNH3 - RB(276)*XH2O*XNH2                             
      W(277) = RF(277)*XO*XNH3 - RB(277)*XOH*XNH2                               

      RETURN
      END

