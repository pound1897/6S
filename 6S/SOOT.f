      subroutine   soot
      real ph,phr
      integer i,j
      common /sixs_aerbas/ ph(10,83)
      dimension phr(10,83)
c
c    model: soot
c
            DATA ((PHR(I,J),J=1,83),I=01,01) /
     *0.4897E+00,0.4896E+00,0.4890E+00,0.4881E+00,0.4867E+00,0.4849E+00,
     *0.4827E+00,0.4802E+00,0.4773E+00,0.4743E+00,0.4709E+00,0.4675E+00,
     *0.4638E+00,0.4601E+00,0.4563E+00,0.4526E+00,0.4489E+00,0.4453E+00,
     *0.4419E+00,0.4388E+00,0.4359E+00,0.4334E+00,0.4312E+00,0.4296E+00,
     *0.4285E+00,0.4281E+00,0.4283E+00,0.4293E+00,0.4312E+00,0.4341E+00,
     *0.4380E+00,0.4430E+00,0.4494E+00,0.4571E+00,0.4663E+00,0.4771E+00,
     *0.4896E+00,0.5041E+00,0.5206E+00,0.5392E+00,0.5603E+00,0.5717E+00,
     *0.5838E+00,0.6101E+00,0.6392E+00,0.6714E+00,0.7069E+00,0.7459E+00,
     *0.7886E+00,0.8352E+00,0.8860E+00,0.9411E+00,0.1001E+01,0.1065E+01,
     *0.1135E+01,0.1210E+01,0.1290E+01,0.1376E+01,0.1468E+01,0.1566E+01,
     *0.1670E+01,0.1781E+01,0.1897E+01,0.2019E+01,0.2148E+01,0.2282E+01,
     *0.2421E+01,0.2565E+01,0.2713E+01,0.2865E+01,0.3019E+01,0.3173E+01,
     *0.3327E+01,0.3479E+01,0.3625E+01,0.3765E+01,0.3894E+01,0.4011E+01,
     *0.4111E+01,0.4192E+01,0.4250E+01,0.4284E+01,0.4292E+01/
            DATA ((PHR(I,J),J=1,83),I=02,02) /
     *0.5620E+00,0.5618E+00,0.5611E+00,0.5599E+00,0.5582E+00,0.5560E+00,
     *0.5533E+00,0.5502E+00,0.5467E+00,0.5428E+00,0.5387E+00,0.5342E+00,
     *0.5295E+00,0.5246E+00,0.5197E+00,0.5146E+00,0.5096E+00,0.5046E+00,
     *0.4998E+00,0.4951E+00,0.4907E+00,0.4866E+00,0.4829E+00,0.4797E+00,
     *0.4771E+00,0.4751E+00,0.4738E+00,0.4734E+00,0.4738E+00,0.4753E+00,
     *0.4779E+00,0.4817E+00,0.4868E+00,0.4934E+00,0.5016E+00,0.5114E+00,
     *0.5231E+00,0.5367E+00,0.5524E+00,0.5704E+00,0.5908E+00,0.6019E+00,
     *0.6137E+00,0.6393E+00,0.6678E+00,0.6993E+00,0.7340E+00,0.7720E+00,
     *0.8136E+00,0.8589E+00,0.9081E+00,0.9613E+00,0.1019E+01,0.1080E+01,
     *0.1147E+01,0.1218E+01,0.1293E+01,0.1373E+01,0.1459E+01,0.1549E+01,
     *0.1643E+01,0.1743E+01,0.1847E+01,0.1956E+01,0.2069E+01,0.2185E+01,
     *0.2305E+01,0.2428E+01,0.2553E+01,0.2679E+01,0.2806E+01,0.2931E+01,
     *0.3055E+01,0.3174E+01,0.3289E+01,0.3396E+01,0.3495E+01,0.3582E+01,
     *0.3656E+01,0.3716E+01,0.3758E+01,0.3782E+01,0.3788E+01/
            DATA ((PHR(I,J),J=1,83),I=03,03) /
     *0.5834E+00,0.5832E+00,0.5825E+00,0.5813E+00,0.5795E+00,0.5771E+00,
     *0.5743E+00,0.5710E+00,0.5673E+00,0.5632E+00,0.5587E+00,0.5540E+00,
     *0.5490E+00,0.5438E+00,0.5384E+00,0.5330E+00,0.5275E+00,0.5221E+00,
     *0.5168E+00,0.5117E+00,0.5068E+00,0.5023E+00,0.4981E+00,0.4944E+00,
     *0.4913E+00,0.4889E+00,0.4871E+00,0.4862E+00,0.4862E+00,0.4872E+00,
     *0.4894E+00,0.4928E+00,0.4975E+00,0.5037E+00,0.5115E+00,0.5210E+00,
     *0.5324E+00,0.5457E+00,0.5611E+00,0.5788E+00,0.5988E+00,0.6098E+00,
     *0.6215E+00,0.6468E+00,0.6749E+00,0.7061E+00,0.7405E+00,0.7781E+00,
     *0.8193E+00,0.8641E+00,0.9127E+00,0.9652E+00,0.1022E+01,0.1083E+01,
     *0.1148E+01,0.1217E+01,0.1291E+01,0.1370E+01,0.1453E+01,0.1541E+01,
     *0.1633E+01,0.1730E+01,0.1831E+01,0.1936E+01,0.2045E+01,0.2157E+01,
     *0.2272E+01,0.2390E+01,0.2509E+01,0.2629E+01,0.2749E+01,0.2867E+01,
     *0.2984E+01,0.3096E+01,0.3203E+01,0.3304E+01,0.3395E+01,0.3476E+01,
     *0.3545E+01,0.3599E+01,0.3638E+01,0.3660E+01,0.3666E+01/
            DATA ((PHR(I,J),J=1,83),I=04,04) /
     *0.6060E+00,0.6059E+00,0.6051E+00,0.6038E+00,0.6019E+00,0.5994E+00,
     *0.5964E+00,0.5929E+00,0.5889E+00,0.5846E+00,0.5798E+00,0.5747E+00,
     *0.5693E+00,0.5637E+00,0.5580E+00,0.5521E+00,0.5462E+00,0.5403E+00,
     *0.5345E+00,0.5289E+00,0.5235E+00,0.5185E+00,0.5138E+00,0.5096E+00,
     *0.5059E+00,0.5029E+00,0.5007E+00,0.4993E+00,0.4988E+00,0.4993E+00,
     *0.5010E+00,0.5040E+00,0.5083E+00,0.5142E+00,0.5216E+00,0.5307E+00,
     *0.5418E+00,0.5548E+00,0.5699E+00,0.5873E+00,0.6071E+00,0.6180E+00,
     *0.6295E+00,0.6546E+00,0.6825E+00,0.7134E+00,0.7474E+00,0.7848E+00,
     *0.8255E+00,0.8699E+00,0.9179E+00,0.9698E+00,0.1026E+01,0.1085E+01,
     *0.1150E+01,0.1218E+01,0.1290E+01,0.1367E+01,0.1448E+01,0.1534E+01,
     *0.1623E+01,0.1717E+01,0.1815E+01,0.1916E+01,0.2020E+01,0.2128E+01,
     *0.2237E+01,0.2349E+01,0.2462E+01,0.2576E+01,0.2688E+01,0.2800E+01,
     *0.2909E+01,0.3013E+01,0.3113E+01,0.3206E+01,0.3290E+01,0.3364E+01,
     *0.3427E+01,0.3477E+01,0.3512E+01,0.3532E+01,0.3537E+01/
            DATA ((PHR(I,J),J=1,83),I=05,05) /
     *0.6604E+00,0.6602E+00,0.6593E+00,0.6578E+00,0.6556E+00,0.6528E+00,
     *0.6494E+00,0.6454E+00,0.6409E+00,0.6358E+00,0.6304E+00,0.6245E+00,
     *0.6182E+00,0.6117E+00,0.6050E+00,0.5981E+00,0.5911E+00,0.5841E+00,
     *0.5771E+00,0.5703E+00,0.5636E+00,0.5573E+00,0.5513E+00,0.5458E+00,
     *0.5409E+00,0.5366E+00,0.5331E+00,0.5305E+00,0.5288E+00,0.5281E+00,
     *0.5287E+00,0.5305E+00,0.5338E+00,0.5385E+00,0.5450E+00,0.5532E+00,
     *0.5633E+00,0.5754E+00,0.5897E+00,0.6062E+00,0.6252E+00,0.6356E+00,
     *0.6467E+00,0.6710E+00,0.6980E+00,0.7280E+00,0.7610E+00,0.7972E+00,
     *0.8367E+00,0.8797E+00,0.9261E+00,0.9762E+00,0.1030E+01,0.1087E+01,
     *0.1149E+01,0.1214E+01,0.1283E+01,0.1355E+01,0.1432E+01,0.1512E+01,
     *0.1595E+01,0.1682E+01,0.1772E+01,0.1865E+01,0.1961E+01,0.2058E+01,
     *0.2157E+01,0.2257E+01,0.2358E+01,0.2458E+01,0.2557E+01,0.2654E+01,
     *0.2748E+01,0.2838E+01,0.2923E+01,0.3001E+01,0.3072E+01,0.3134E+01,
     *0.3187E+01,0.3228E+01,0.3257E+01,0.3273E+01,0.3277E+01/
            DATA ((PHR(I,J),J=1,83),I=06,06) /
     *0.6993E+00,0.6991E+00,0.6982E+00,0.6965E+00,0.6942E+00,0.6911E+00,
     *0.6874E+00,0.6830E+00,0.6781E+00,0.6726E+00,0.6666E+00,0.6601E+00,
     *0.6533E+00,0.6461E+00,0.6387E+00,0.6310E+00,0.6232E+00,0.6154E+00,
     *0.6076E+00,0.5998E+00,0.5923E+00,0.5851E+00,0.5782E+00,0.5717E+00,
     *0.5659E+00,0.5607E+00,0.5562E+00,0.5526E+00,0.5500E+00,0.5485E+00,
     *0.5482E+00,0.5491E+00,0.5515E+00,0.5555E+00,0.5611E+00,0.5686E+00,
     *0.5779E+00,0.5893E+00,0.6028E+00,0.6187E+00,0.6369E+00,0.6470E+00,
     *0.6577E+00,0.6812E+00,0.7074E+00,0.7366E+00,0.7687E+00,0.8040E+00,
     *0.8425E+00,0.8843E+00,0.9295E+00,0.9781E+00,0.1030E+01,0.1086E+01,
     *0.1145E+01,0.1208E+01,0.1274E+01,0.1344E+01,0.1417E+01,0.1494E+01,
     *0.1573E+01,0.1656E+01,0.1741E+01,0.1828E+01,0.1918E+01,0.2009E+01,
     *0.2101E+01,0.2194E+01,0.2287E+01,0.2380E+01,0.2470E+01,0.2559E+01,
     *0.2645E+01,0.2726E+01,0.2803E+01,0.2873E+01,0.2937E+01,0.2992E+01,
     *0.3038E+01,0.3075E+01,0.3100E+01,0.3115E+01,0.3118E+01/
            DATA ((PHR(I,J),J=1,83),I=07,07) /
     *0.7916E+00,0.7914E+00,0.7903E+00,0.7883E+00,0.7855E+00,0.7818E+00,
     *0.7773E+00,0.7721E+00,0.7662E+00,0.7595E+00,0.7522E+00,0.7444E+00,
     *0.7360E+00,0.7272E+00,0.7180E+00,0.7085E+00,0.6988E+00,0.6889E+00,
     *0.6790E+00,0.6692E+00,0.6595E+00,0.6500E+00,0.6408E+00,0.6321E+00,
     *0.6239E+00,0.6164E+00,0.6097E+00,0.6038E+00,0.5989E+00,0.5952E+00,
     *0.5926E+00,0.5915E+00,0.5918E+00,0.5936E+00,0.5972E+00,0.6027E+00,
     *0.6101E+00,0.6195E+00,0.6311E+00,0.6451E+00,0.6614E+00,0.6705E+00,
     *0.6803E+00,0.7017E+00,0.7259E+00,0.7529E+00,0.7828E+00,0.8156E+00,
     *0.8514E+00,0.8903E+00,0.9323E+00,0.9774E+00,0.1026E+01,0.1077E+01,
     *0.1131E+01,0.1189E+01,0.1249E+01,0.1312E+01,0.1378E+01,0.1447E+01,
     *0.1518E+01,0.1590E+01,0.1665E+01,0.1741E+01,0.1819E+01,0.1897E+01,
     *0.1976E+01,0.2054E+01,0.2132E+01,0.2209E+01,0.2284E+01,0.2356E+01,
     *0.2426E+01,0.2491E+01,0.2552E+01,0.2607E+01,0.2657E+01,0.2700E+01,
     *0.2736E+01,0.2764E+01,0.2783E+01,0.2795E+01,0.2797E+01/
            DATA ((PHR(I,J),J=1,83),I=08,08) /
     *0.1041E+01,0.1040E+01,0.1038E+01,0.1036E+01,0.1031E+01,0.1026E+01,
     *0.1019E+01,0.1011E+01,0.1002E+01,0.9924E+00,0.9814E+00,0.9694E+00,
     *0.9566E+00,0.9431E+00,0.9288E+00,0.9140E+00,0.8988E+00,0.8832E+00,
     *0.8673E+00,0.8513E+00,0.8353E+00,0.8194E+00,0.8038E+00,0.7885E+00,
     *0.7737E+00,0.7596E+00,0.7462E+00,0.7338E+00,0.7223E+00,0.7121E+00,
     *0.7031E+00,0.6955E+00,0.6895E+00,0.6852E+00,0.6827E+00,0.6820E+00,
     *0.6833E+00,0.6868E+00,0.6924E+00,0.7003E+00,0.7105E+00,0.7165E+00,
     *0.7232E+00,0.7383E+00,0.7559E+00,0.7760E+00,0.7987E+00,0.8240E+00,
     *0.8518E+00,0.8821E+00,0.9149E+00,0.9501E+00,0.9877E+00,0.1028E+01,
     *0.1069E+01,0.1113E+01,0.1159E+01,0.1207E+01,0.1256E+01,0.1306E+01,
     *0.1358E+01,0.1410E+01,0.1463E+01,0.1517E+01,0.1570E+01,0.1623E+01,
     *0.1676E+01,0.1727E+01,0.1778E+01,0.1827E+01,0.1873E+01,0.1918E+01,
     *0.1960E+01,0.1999E+01,0.2035E+01,0.2067E+01,0.2096E+01,0.2120E+01,
     *0.2140E+01,0.2156E+01,0.2167E+01,0.2173E+01,0.2174E+01/
            DATA ((PHR(I,J),J=1,83),I=09,09) /
     *0.1182E+01,0.1181E+01,0.1179E+01,0.1176E+01,0.1171E+01,0.1164E+01,
     *0.1156E+01,0.1147E+01,0.1136E+01,0.1124E+01,0.1110E+01,0.1096E+01,
     *0.1080E+01,0.1064E+01,0.1046E+01,0.1028E+01,0.1009E+01,0.9903E+00,
     *0.9708E+00,0.9510E+00,0.9312E+00,0.9114E+00,0.8919E+00,0.8726E+00,
     *0.8539E+00,0.8357E+00,0.8184E+00,0.8019E+00,0.7866E+00,0.7724E+00,
     *0.7595E+00,0.7481E+00,0.7383E+00,0.7302E+00,0.7239E+00,0.7195E+00,
     *0.7171E+00,0.7168E+00,0.7188E+00,0.7229E+00,0.7294E+00,0.7335E+00,
     *0.7382E+00,0.7494E+00,0.7630E+00,0.7790E+00,0.7974E+00,0.8182E+00,
     *0.8414E+00,0.8668E+00,0.8944E+00,0.9242E+00,0.9561E+00,0.9898E+00,
     *0.1025E+01,0.1063E+01,0.1101E+01,0.1141E+01,0.1183E+01,0.1225E+01,
     *0.1268E+01,0.1311E+01,0.1355E+01,0.1399E+01,0.1442E+01,0.1485E+01,
     *0.1528E+01,0.1569E+01,0.1609E+01,0.1648E+01,0.1685E+01,0.1720E+01,
     *0.1753E+01,0.1783E+01,0.1811E+01,0.1836E+01,0.1858E+01,0.1876E+01,
     *0.1891E+01,0.1903E+01,0.1911E+01,0.1916E+01,0.1917E+01/
            DATA ((PHR(I,J),J=1,83),I=10,10) /
     *0.1325E+01,0.1324E+01,0.1322E+01,0.1318E+01,0.1312E+01,0.1304E+01,
     *0.1294E+01,0.1283E+01,0.1270E+01,0.1256E+01,0.1240E+01,0.1222E+01,
     *0.1204E+01,0.1184E+01,0.1163E+01,0.1142E+01,0.1119E+01,0.1096E+01,
     *0.1073E+01,0.1049E+01,0.1025E+01,0.1001E+01,0.9776E+00,0.9541E+00,
     *0.9312E+00,0.9088E+00,0.8872E+00,0.8666E+00,0.8471E+00,0.8287E+00,
     *0.8118E+00,0.7963E+00,0.7825E+00,0.7704E+00,0.7602E+00,0.7519E+00,
     *0.7457E+00,0.7415E+00,0.7396E+00,0.7399E+00,0.7424E+00,0.7446E+00,
     *0.7473E+00,0.7545E+00,0.7640E+00,0.7758E+00,0.7899E+00,0.8063E+00,
     *0.8248E+00,0.8455E+00,0.8681E+00,0.8928E+00,0.9192E+00,0.9473E+00,
     *0.9771E+00,0.1008E+01,0.1041E+01,0.1074E+01,0.1109E+01,0.1144E+01,
     *0.1179E+01,0.1215E+01,0.1252E+01,0.1288E+01,0.1324E+01,0.1359E+01,
     *0.1393E+01,0.1427E+01,0.1460E+01,0.1491E+01,0.1521E+01,0.1549E+01,
     *0.1575E+01,0.1599E+01,0.1622E+01,0.1641E+01,0.1658E+01,0.1673E+01,
     *0.1685E+01,0.1694E+01,0.1701E+01,0.1704E+01,0.1705E+01/
      do 1 i=1,10
      do 1 j=1,83
      ph(i,j)=phr(i,j)
    1 continue
      return
      end
