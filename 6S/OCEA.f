      subroutine   ocea
      common /sixs_aerbas/ ph(10,83)
      real phr(10,83),ph
      integer i,j
c
c    model: oceanic
c
            DATA ((PHR(I,J),J=1,83),I=01,01) /
     *0.7855E+00,0.6283E+00,0.5465E+00,0.4693E+00,0.4153E+00,0.3917E+00,
     *0.3657E+00,0.3378E+00,0.3161E+00,0.3025E+00,0.2972E+00,0.2990E+00,
     *0.3055E+00,0.3118E+00,0.3059E+00,0.2715E+00,0.2118E+00,0.1585E+00,
     *0.1230E+00,0.9913E-01,0.8327E-01,0.7292E-01,0.6585E-01,0.6171E-01,
     *0.5883E-01,0.5780E-01,0.5791E-01,0.5893E-01,0.6144E-01,0.6406E-01,
     *0.6717E-01,0.6966E-01,0.7130E-01,0.7291E-01,0.7434E-01,0.7626E-01,
     *0.7847E-01,0.8190E-01,0.8583E-01,0.9044E-01,0.9709E-01,0.1006E+00,
     *0.1045E+00,0.1128E+00,0.1239E+00,0.1360E+00,0.1497E+00,0.1667E+00,
     *0.1856E+00,0.2070E+00,0.2323E+00,0.2615E+00,0.2948E+00,0.3326E+00,
     *0.3772E+00,0.4263E+00,0.4840E+00,0.5492E+00,0.6242E+00,0.7103E+00,
     *0.8075E+00,0.9192E+00,0.1046E+01,0.1190E+01,0.1354E+01,0.1541E+01,
     *0.1756E+01,0.2002E+01,0.2277E+01,0.2603E+01,0.2976E+01,0.3416E+01,
     *0.3931E+01,0.4563E+01,0.5372E+01,0.6490E+01,0.8191E+01,0.1111E+02,
     *0.1692E+02,0.3097E+02,0.7524E+02,0.2992E+03,0.1697E+04/
            DATA ((PHR(I,J),J=1,83),I=02,02) /
     *0.7129E+00,0.5739E+00,0.5059E+00,0.4429E+00,0.4035E+00,0.3898E+00,
     *0.3678E+00,0.3416E+00,0.3195E+00,0.3042E+00,0.2975E+00,0.2961E+00,
     *0.2987E+00,0.2994E+00,0.2909E+00,0.2614E+00,0.2134E+00,0.1670E+00,
     *0.1336E+00,0.1100E+00,0.9363E-01,0.8252E-01,0.7480E-01,0.6967E-01,
     *0.6621E-01,0.6499E-01,0.6438E-01,0.6506E-01,0.6656E-01,0.6880E-01,
     *0.7108E-01,0.7332E-01,0.7497E-01,0.7681E-01,0.7860E-01,0.8093E-01,
     *0.8357E-01,0.8723E-01,0.9184E-01,0.9665E-01,0.1036E+00,0.1075E+00,
     *0.1112E+00,0.1200E+00,0.1316E+00,0.1436E+00,0.1580E+00,0.1748E+00,
     *0.1937E+00,0.2154E+00,0.2413E+00,0.2704E+00,0.3031E+00,0.3421E+00,
     *0.3856E+00,0.4356E+00,0.4928E+00,0.5586E+00,0.6333E+00,0.7196E+00,
     *0.8188E+00,0.9313E+00,0.1060E+01,0.1208E+01,0.1375E+01,0.1568E+01,
     *0.1791E+01,0.2047E+01,0.2340E+01,0.2679E+01,0.3075E+01,0.3547E+01,
     *0.4107E+01,0.4805E+01,0.5714E+01,0.6981E+01,0.8889E+01,0.1212E+02,
     *0.1839E+02,0.3283E+02,0.7515E+02,0.2626E+03,0.1134E+04/
            DATA ((PHR(I,J),J=1,83),I=03,03) /
     *0.6966E+00,0.5607E+00,0.4902E+00,0.4336E+00,0.3978E+00,0.3866E+00,
     *0.3674E+00,0.3412E+00,0.3187E+00,0.3039E+00,0.2960E+00,0.2945E+00,
     *0.2960E+00,0.2961E+00,0.2874E+00,0.2591E+00,0.2133E+00,0.1692E+00,
     *0.1362E+00,0.1129E+00,0.9630E-01,0.8484E-01,0.7707E-01,0.7190E-01,
     *0.6854E-01,0.6653E-01,0.6597E-01,0.6668E-01,0.6812E-01,0.7009E-01,
     *0.7216E-01,0.7425E-01,0.7580E-01,0.7758E-01,0.7959E-01,0.8174E-01,
     *0.8490E-01,0.8852E-01,0.9294E-01,0.9864E-01,0.1048E+00,0.1084E+00,
     *0.1128E+00,0.1220E+00,0.1325E+00,0.1453E+00,0.1596E+00,0.1762E+00,
     *0.1959E+00,0.2177E+00,0.2428E+00,0.2725E+00,0.3055E+00,0.3440E+00,
     *0.3882E+00,0.4382E+00,0.4953E+00,0.5613E+00,0.6365E+00,0.7225E+00,
     *0.8218E+00,0.9344E+00,0.1065E+01,0.1212E+01,0.1381E+01,0.1577E+01,
     *0.1801E+01,0.2059E+01,0.2360E+01,0.2701E+01,0.3107E+01,0.3586E+01,
     *0.4166E+01,0.4885E+01,0.5821E+01,0.7115E+01,0.9088E+01,0.1241E+02,
     *0.1877E+02,0.3323E+02,0.7480E+02,0.2523E+03,0.1018E+04/
            DATA ((PHR(I,J),J=1,83),I=04,04) /
     *0.6774E+00,0.5476E+00,0.4775E+00,0.4252E+00,0.3937E+00,0.3855E+00,
     *0.3684E+00,0.3432E+00,0.3209E+00,0.3059E+00,0.2974E+00,0.2950E+00,
     *0.2951E+00,0.2935E+00,0.2832E+00,0.2550E+00,0.2114E+00,0.1697E+00,
     *0.1380E+00,0.1153E+00,0.9882E-01,0.8737E-01,0.7952E-01,0.7423E-01,
     *0.7074E-01,0.6859E-01,0.6788E-01,0.6842E-01,0.6969E-01,0.7150E-01,
     *0.7349E-01,0.7557E-01,0.7720E-01,0.7911E-01,0.8125E-01,0.8356E-01,
     *0.8685E-01,0.9062E-01,0.9516E-01,0.1010E+00,0.1073E+00,0.1109E+00,
     *0.1154E+00,0.1247E+00,0.1352E+00,0.1482E+00,0.1626E+00,0.1793E+00,
     *0.1991E+00,0.2210E+00,0.2462E+00,0.2760E+00,0.3091E+00,0.3477E+00,
     *0.3920E+00,0.4422E+00,0.4994E+00,0.5656E+00,0.6410E+00,0.7275E+00,
     *0.8272E+00,0.9405E+00,0.1071E+01,0.1220E+01,0.1391E+01,0.1588E+01,
     *0.1815E+01,0.2077E+01,0.2382E+01,0.2731E+01,0.3145E+01,0.3636E+01,
     *0.4233E+01,0.4974E+01,0.5942E+01,0.7282E+01,0.9319E+01,0.1273E+02,
     *0.1919E+02,0.3364E+02,0.7414E+02,0.2397E+03,0.8914E+03/
            DATA ((PHR(I,J),J=1,83),I=05,05) /
     *0.6153E+00,0.5058E+00,0.4382E+00,0.3950E+00,0.3738E+00,0.3731E+00,
     *0.3585E+00,0.3354E+00,0.3139E+00,0.2983E+00,0.2892E+00,0.2849E+00,
     *0.2832E+00,0.2800E+00,0.2703E+00,0.2469E+00,0.2112E+00,0.1741E+00,
     *0.1442E+00,0.1219E+00,0.1054E+00,0.9356E-01,0.8531E-01,0.7966E-01,
     *0.7561E-01,0.7323E-01,0.7198E-01,0.7214E-01,0.7291E-01,0.7415E-01,
     *0.7601E-01,0.7747E-01,0.7901E-01,0.8091E-01,0.8293E-01,0.8564E-01,
     *0.8906E-01,0.9289E-01,0.9788E-01,0.1033E+00,0.1102E+00,0.1141E+00,
     *0.1181E+00,0.1275E+00,0.1385E+00,0.1511E+00,0.1660E+00,0.1823E+00,
     *0.2018E+00,0.2241E+00,0.2491E+00,0.2784E+00,0.3123E+00,0.3503E+00,
     *0.3942E+00,0.4451E+00,0.5020E+00,0.5684E+00,0.6448E+00,0.7319E+00,
     *0.8325E+00,0.9481E+00,0.1081E+01,0.1234E+01,0.1409E+01,0.1612E+01,
     *0.1846E+01,0.2118E+01,0.2440E+01,0.2809E+01,0.3249E+01,0.3773E+01,
     *0.4413E+01,0.5211E+01,0.6259E+01,0.7710E+01,0.9888E+01,0.1347E+02,
     *0.2009E+02,0.3435E+02,0.7217E+02,0.2130E+03,0.6728E+03/
            DATA ((PHR(I,J),J=1,83),I=06,06) /
     *0.5916E+00,0.4877E+00,0.4171E+00,0.3786E+00,0.3632E+00,0.3654E+00,
     *0.3546E+00,0.3335E+00,0.3124E+00,0.2967E+00,0.2869E+00,0.2822E+00,
     *0.2792E+00,0.2744E+00,0.2635E+00,0.2413E+00,0.2085E+00,0.1740E+00,
     *0.1459E+00,0.1244E+00,0.1084E+00,0.9682E-01,0.8822E-01,0.8243E-01,
     *0.7835E-01,0.7606E-01,0.7463E-01,0.7441E-01,0.7473E-01,0.7609E-01,
     *0.7739E-01,0.7905E-01,0.8078E-01,0.8256E-01,0.8474E-01,0.8745E-01,
     *0.9082E-01,0.9490E-01,0.9996E-01,0.1057E+00,0.1127E+00,0.1166E+00,
     *0.1207E+00,0.1301E+00,0.1412E+00,0.1539E+00,0.1686E+00,0.1858E+00,
     *0.2048E+00,0.2270E+00,0.2528E+00,0.2818E+00,0.3154E+00,0.3545E+00,
     *0.3980E+00,0.4487E+00,0.5067E+00,0.5728E+00,0.6491E+00,0.7374E+00,
     *0.8386E+00,0.9547E+00,0.1090E+01,0.1244E+01,0.1423E+01,0.1630E+01,
     *0.1870E+01,0.2149E+01,0.2477E+01,0.2862E+01,0.3316E+01,0.3862E+01,
     *0.4527E+01,0.5365E+01,0.6458E+01,0.7974E+01,0.1023E+02,0.1390E+02,
     *0.2058E+02,0.3459E+02,0.7042E+02,0.1961E+03,0.5608E+03/
            DATA ((PHR(I,J),J=1,83),I=07,07) /
     *0.5164E+00,0.4330E+00,0.3650E+00,0.3341E+00,0.3313E+00,0.3413E+00,
     *0.3356E+00,0.3182E+00,0.2998E+00,0.2844E+00,0.2744E+00,0.2677E+00,
     *0.2626E+00,0.2560E+00,0.2453E+00,0.2267E+00,0.2009E+00,0.1730E+00,
     *0.1485E+00,0.1291E+00,0.1141E+00,0.1028E+00,0.9425E-01,0.8828E-01,
     *0.8375E-01,0.8105E-01,0.7927E-01,0.7843E-01,0.7860E-01,0.7925E-01,
     *0.8010E-01,0.8165E-01,0.8331E-01,0.8499E-01,0.8754E-01,0.9034E-01,
     *0.9390E-01,0.9825E-01,0.1034E+00,0.1093E+00,0.1164E+00,0.1203E+00,
     *0.1246E+00,0.1342E+00,0.1452E+00,0.1582E+00,0.1728E+00,0.1896E+00,
     *0.2094E+00,0.2310E+00,0.2569E+00,0.2863E+00,0.3195E+00,0.3587E+00,
     *0.4030E+00,0.4534E+00,0.5122E+00,0.5794E+00,0.6565E+00,0.7463E+00,
     *0.8505E+00,0.9697E+00,0.1109E+01,0.1270E+01,0.1457E+01,0.1674E+01,
     *0.1929E+01,0.2226E+01,0.2578E+01,0.2997E+01,0.3495E+01,0.4096E+01,
     *0.4831E+01,0.5758E+01,0.6967E+01,0.8629E+01,0.1105E+02,0.1487E+02,
     *0.2152E+02,0.3465E+02,0.6548E+02,0.1595E+03,0.3700E+03/
            DATA ((PHR(I,J),J=1,83),I=08,08) /
     *0.3257E+00,0.2888E+00,0.2378E+00,0.2215E+00,0.2345E+00,0.2532E+00,
     *0.2578E+00,0.2504E+00,0.2390E+00,0.2282E+00,0.2194E+00,0.2123E+00,
     *0.2059E+00,0.1991E+00,0.1906E+00,0.1797E+00,0.1665E+00,0.1520E+00,
     *0.1379E+00,0.1254E+00,0.1147E+00,0.1061E+00,0.9917E-01,0.9373E-01,
     *0.8960E-01,0.8656E-01,0.8438E-01,0.8306E-01,0.8243E-01,0.8240E-01,
     *0.8294E-01,0.8394E-01,0.8543E-01,0.8740E-01,0.8990E-01,0.9302E-01,
     *0.9681E-01,0.1013E+00,0.1067E+00,0.1129E+00,0.1200E+00,0.1240E+00,
     *0.1283E+00,0.1379E+00,0.1490E+00,0.1618E+00,0.1764E+00,0.1932E+00,
     *0.2124E+00,0.2345E+00,0.2599E+00,0.2892E+00,0.3231E+00,0.3622E+00,
     *0.4072E+00,0.4593E+00,0.5195E+00,0.5895E+00,0.6711E+00,0.7664E+00,
     *0.8781E+00,0.1009E+01,0.1163E+01,0.1343E+01,0.1556E+01,0.1808E+01,
     *0.2107E+01,0.2464E+01,0.2891E+01,0.3405E+01,0.4025E+01,0.4779E+01,
     *0.5707E+01,0.6863E+01,0.8338E+01,0.1027E+02,0.1291E+02,0.1670E+02,
     *0.2248E+02,0.3211E+02,0.5001E+02,0.8772E+02,0.1334E+03/
            DATA ((PHR(I,J),J=1,83),I=09,09) /
     *0.2139E+00,0.1949E+00,0.1618E+00,0.1541E+00,0.1685E+00,0.1828E+00,
     *0.1856E+00,0.1800E+00,0.1718E+00,0.1642E+00,0.1581E+00,0.1534E+00,
     *0.1495E+00,0.1460E+00,0.1421E+00,0.1375E+00,0.1318E+00,0.1252E+00,
     *0.1178E+00,0.1105E+00,0.1036E+00,0.9754E-01,0.9237E-01,0.8811E-01,
     *0.8468E-01,0.8198E-01,0.7994E-01,0.7852E-01,0.7768E-01,0.7741E-01,
     *0.7767E-01,0.7843E-01,0.7969E-01,0.8144E-01,0.8373E-01,0.8662E-01,
     *0.9014E-01,0.9438E-01,0.9939E-01,0.1052E+00,0.1120E+00,0.1158E+00,
     *0.1198E+00,0.1289E+00,0.1394E+00,0.1514E+00,0.1653E+00,0.1813E+00,
     *0.1997E+00,0.2208E+00,0.2453E+00,0.2736E+00,0.3064E+00,0.3444E+00,
     *0.3886E+00,0.4400E+00,0.5000E+00,0.5703E+00,0.6528E+00,0.7502E+00,
     *0.8652E+00,0.1001E+01,0.1163E+01,0.1355E+01,0.1584E+01,0.1859E+01,
     *0.2188E+01,0.2586E+01,0.3067E+01,0.3649E+01,0.4358E+01,0.5222E+01,
     *0.6282E+01,0.7594E+01,0.9235E+01,0.1132E+02,0.1404E+02,0.1768E+02,
     *0.2278E+02,0.3033E+02,0.4233E+02,0.6237E+02,0.7953E+02/
            DATA ((PHR(I,J),J=1,83),I=10,10) /
     *0.2110E+00,0.2025E+00,0.1832E+00,0.1730E+00,0.1773E+00,0.1912E+00,
     *0.2055E+00,0.2138E+00,0.2152E+00,0.2113E+00,0.2040E+00,0.1946E+00,
     *0.1842E+00,0.1734E+00,0.1627E+00,0.1524E+00,0.1429E+00,0.1344E+00,
     *0.1268E+00,0.1203E+00,0.1149E+00,0.1104E+00,0.1068E+00,0.1040E+00,
     *0.1019E+00,0.1006E+00,0.9982E-01,0.9972E-01,0.1003E+00,0.1014E+00,
     *0.1031E+00,0.1054E+00,0.1084E+00,0.1119E+00,0.1162E+00,0.1212E+00,
     *0.1271E+00,0.1338E+00,0.1415E+00,0.1503E+00,0.1603E+00,0.1658E+00,
     *0.1717E+00,0.1847E+00,0.1995E+00,0.2163E+00,0.2354E+00,0.2571E+00,
     *0.2818E+00,0.3100E+00,0.3422E+00,0.3792E+00,0.4216E+00,0.4702E+00,
     *0.5261E+00,0.5903E+00,0.6644E+00,0.7500E+00,0.8493E+00,0.9645E+00,
     *0.1098E+01,0.1254E+01,0.1436E+01,0.1649E+01,0.1897E+01,0.2189E+01,
     *0.2531E+01,0.2934E+01,0.3408E+01,0.3968E+01,0.4630E+01,0.5415E+01,
     *0.6348E+01,0.7463E+01,0.8805E+01,0.1044E+02,0.1244E+02,0.1495E+02,
     *0.1816E+02,0.2237E+02,0.2799E+02,0.3517E+02,0.3934E+02/
      do 1 i=1,10
      do 1 j=1,83
      ph(i,j)=phr(i,j)
    1 continue
      return
      end
