      subroutine   SuperView1(iwa)
      real s,wlinf,wlsup
      common /sixs_ffu/ s(1501),wlinf,wlsup
      real sr(5,1501),wli(5),wls(5)
      integer iwa,l,i

c    pan spectral band of Superview1-01

      data (sr(1,l),l=1,1501)/  60*0.,
     a 0.0041,0.0043,0.0046,0.0049,0.0054,0.0057,0.0063,
     a 0.0065,0.0069,0.0072,0.0072,0.0071,0.0072,0.0075,
     a 0.0077,0.0078,0.0079,0.0101,0.0183,0.0476,0.0904,
     a 0.1140,0.1244,0.1328,0.1422,0.1510,0.1590,0.1673,
     a 0.1745,0.1799,0.1834,0.1887,0.1977,0.2061,0.2152,
     a 0.2250,0.2353,0.2454,0.2522,0.2577,0.2611,0.2656,
     a 0.2725,0.2829,0.2953,0.3103,0.3242,0.3364,0.3487,
     a 0.3561,0.3648,0.3728,0.3813,0.3894,0.3962,0.4029,
     a 0.4097,0.4137,0.4154,0.4170,0.4176,0.4175,0.4164,
     a 0.4145,0.4120,0.4123,0.4160,0.4214,0.4267,0.4329,
     a 0.4380,0.4411,0.4394,0.4341,0.4278,0.4249,0.4235,
     a 0.4238,0.4255,0.4299,0.4350,0.4440,0.4540,0.4669,
     a 0.4832,0.5043,0.5273,0.5518,0.5774,0.5991,0.6181,
     a 0.6316,0.6376,0.6439,0.6521,0.6636,0.6796,0.6969,
     a 0.7153,0.7307,0.7440,0.7534,0.7568,0.7565,0.7532,
     a 0.7423,0.7275,0.7115,0.6961,0.6865,0.6816,0.6795,
     a 0.6793,0.6834,0.6894,0.6938,0.6951,0.6899,0.6839,
     a 0.6783,0.6725,0.6646,0.6574,0.6538,0.6563,0.6638,
     a 0.6742,0.6855,0.6913,0.6932,0.6925,0.6907,0.6887,
     a 0.6845,0.6771,0.6693,0.6629,0.6591,0.6579,0.6576,
     a 0.6599,0.6634,0.6712,0.6842,0.6983,0.7163,0.7351,
     a 0.7533,0.7690,0.7794,0.7885,0.7999,0.8136,0.8275,
     a 0.8428,0.8560,0.8691,0.8819,0.8920,0.9005,0.9067,
     a 0.9105,0.9171,0.9260,0.9369,0.9502,0.9620,0.9737,
     a 0.9848,0.9944,0.9997,0.9993,0.9941,0.9893,0.9829,
     a 0.9816,0.9833,0.9823,0.9867,0.9879,0.9873,0.9847,
     a 0.9764,0.9663,0.9579,0.9519,0.9478,0.9475,0.9480,
     a 0.9492,0.9443,0.9345,0.9226,0.9128,0.9113,0.9142,
     a 0.8970,0.8082,0.6351,0.4253,0.2516,0.1440,0.0843,
     a 0.0505,0.0323,0.0218,0.0161,0.0118,0.0101,0.0078,
     a 0.0071,0.0063,0.0061,0.0059,0.0051,0.0051,0.0051,
     a 0.0051,0.0051,0.0051,0.0051,0.0049,0.0041,0.0041,
     a 0.0041,0.0041,0.0041,0.0041,0.0041,0.0041,0.0041,
     a 0.0040,0.0051,0.0051,0.0051,0.0051,0.0051,0.0051,
     a1203*0./
c 
c    1st spectral band of Superview1-01
 
      data (sr(2,l),l=1,1501)/  60*0.,
     a 0.0001,0.0003,0.0002,0.0003,0.0001,0.0001,0.0001,
     a 0.0002,0.0006,0.0012,0.0015,0.0016,0.0018,0.0018,
     a 0.0019,0.0024,0.0039,0.0072,0.0211,0.0816,0.2077,
     a 0.3059,0.3467,0.3739,0.4043,0.4334,0.4612,0.4935,
     a 0.5300,0.5572,0.5716,0.5855,0.6131,0.6477,0.6772,
     a 0.7025,0.7310,0.7656,0.7975,0.8166,0.8224,0.8299,
     a 0.8587,0.9055,0.9457,0.9800,0.9960,0.8956,0.5979,
     a 0.2852,0.1212,0.0549,0.0293,0.0200,0.0148,0.0125,
     a 0.0112,0.0102,0.0108,0.0121,0.0133,0.0159,0.0193,
     a 0.0221,0.0227,0.0223,0.0201,0.0187,0.0185,0.0188,
     a 0.0201,0.0215,0.0217,0.0215,0.0201,0.0186,0.0165,
     a 0.0139,0.0095,0.0054,0.0030,0.0017,0.0011,0.0007,
     a 0.0004,0.0003,0.0002,0.0001,0.0001,0.0001,0.0001,
     a1350*0./
c
c    2nd spectral band of Superview1-01
      data (sr(3,l),l=1,1501)/ 60*0.,
     a 0.0005,0.0008,0.0007,0.0008,0.0002,0.0001,0.0002,
     a 0.0004,0.0010,0.0019,0.0025,0.0022,0.0023,0.0023,
     a 0.0018,0.0019,0.0019,0.0019,0.0023,0.0026,0.0034,
     a 0.0042,0.0042,0.0043,0.0045,0.0048,0.0055,0.0063,
     a 0.0074,0.0076,0.0077,0.0077,0.0087,0.0098,0.0106,
     a 0.0093,0.0081,0.0063,0.0061,0.0056,0.0057,0.0069,
     a 0.0081,0.0108,0.0165,0.0341,0.0892,0.2291,0.4708,
     a 0.6737,0.7498,0.7774,0.8077,0.8427,0.8724,0.9042,
     a 0.9389,0.9723,0.9941,0.9997,0.9914,0.9774,0.9702,
     a 0.9671,0.9642,0.9588,0.9558,0.9558,0.9566,0.9498,
     a 0.9406,0.9330,0.9255,0.9194,0.9116,0.8680,0.6901,
     a 0.3753,0.1276,0.0305,0.0068,0.0019,0.0006,0.0002,
     a1357*0./
c
c    3nd spectral band of Superview1-01
 
      data (sr(4,l),l=1,1501)/ 60*0.,
     a 0.0012,0.0014,0.0016,0.0015,0.0016,0.0018,0.0018,
     a 0.0018,0.0019,0.0020,0.0019,0.0017,0.0017,0.0018,
     a 0.0018,0.0014,0.0013,0.0013,0.0016,0.0030,0.0034,
     a 0.0035,0.0039,0.0049,0.0052,0.0052,0.0047,0.0048,
     a 0.0058,0.0069,0.0069,0.0054,0.0044,0.0045,0.0050,
     a 0.0061,0.0083,0.0082,0.0066,0.0049,0.0047,0.0050,
     a 0.0063,0.0083,0.0097,0.0096,0.0083,0.0079,0.0079,
     a 0.0091,0.0104,0.0117,0.0109,0.0097,0.0090,0.0090,
     a 0.0090,0.0102,0.0115,0.0120,0.0120,0.0118,0.0105,
     a 0.0093,0.0080,0.0070,0.0070,0.0070,0.0080,0.0095,
     a 0.0120,0.0146,0.0151,0.0125,0.0095,0.0073,0.0060,
     a 0.0048,0.0040,0.0040,0.0040,0.0042,0.0050,0.0058,
     a 0.0070,0.0095,0.0130,0.0209,0.0350,0.0715,0.1527,
     a 0.3121,0.5517,0.7635,0.8609,0.8928,0.9173,0.9448,
     a 0.9679,0.9796,0.9844,0.9916,0.9976,0.9994,0.9942,
     a 0.9741,0.9530,0.9377,0.9296,0.9242,0.9141,0.8991,
     a 0.8862,0.8826,0.8615,0.7660,0.5629,0.3357,0.1775,
     a 0.0908,0.0480,0.0282,0.0180,0.0119,0.0090,0.0067,
     a 0.0055,0.0042,0.0040,0.0030,0.0030,0.0030,0.0030,
     a 0.0030,0.0030,0.0030,0.0030,0.0030,0.0035,0.0039,
     a 0.0039,0.0052,0.0059,0.0066,0.0079,0.0089,0.0093,
     a 0.0097,0.0107,0.0107,0.0112,0.0116,0.0127,0.0138,
     a 0.0150,0.0170,0.0194,0.0220,0.0252,0.0284,0.0300,
     a 0.0296,0.0271,0.0246,0.0222,0.0198,0.0174,0.0157,
     a 0.0145,0.0143,0.0135,0.0135,0.0144,0.0145,0.0149,
     a 0.0162,0.0174,0.0176,0.0184,0.0191,0.0193,0.0193,
     a 0.0193,0.0193,0.0183,0.0184,0.0184,0.0176,0.0184,
     a 0.0185,0.0185,0.0192,0.0194,0.0195,0.0199,0.0205,
     a 0.0195,0.0194,0.0190,0.0178,0.0166,0.0151,0.0131,
     a 0.0111,0.0097,0.0075,0.0063,0.0051,0.0039,0.0029,
     a 0.0024,0.0020,0.0019,0.0009,0.0010,0.0010,0.0010,
     a 0.0010,0.0010,0.0010,0.0010,0.0010,0.0010,0.0010,
     a1217*0./
c
c    4th spectral band of Superview1-01
     
      data (sr(5,l),l=1,1501)/ 60*0.,
     a 0.0013,0.0014,0.0014,0.0015,0.0014,0.0015,0.0017,
     a 0.0018,0.0020,0.0021,0.0022,0.0023,0.0024,0.0023,
     a 0.0019,0.0019,0.0020,0.0022,0.0027,0.0028,0.0035,
     a 0.0036,0.0037,0.0043,0.0046,0.0050,0.0057,0.0065,
     a 0.0068,0.0059,0.0058,0.0063,0.0055,0.0044,0.0037,
     a 0.0045,0.0057,0.0077,0.0073,0.0068,0.0079,0.0094,
     a 0.0110,0.0103,0.0091,0.0088,0.0076,0.0071,0.0072,
     a 0.0074,0.0082,0.0090,0.0103,0.0116,0.0124,0.0132,
     a 0.0124,0.0124,0.0124,0.0132,0.0135,0.0119,0.0099,
     a 0.0094,0.0104,0.0127,0.0135,0.0127,0.0094,0.0073,
     a 0.0073,0.0073,0.0083,0.0099,0.0120,0.0125,0.0125,
     a 0.0120,0.0099,0.0086,0.0073,0.0065,0.0073,0.0096,
     a 0.0135,0.0156,0.0140,0.0101,0.0073,0.0073,0.0078,
     a 0.0099,0.0156,0.0236,0.0244,0.0189,0.0135,0.0114,
     a 0.0124,0.0158,0.0228,0.0257,0.0228,0.0168,0.0124,
     a 0.0117,0.0145,0.0213,0.0291,0.0249,0.0151,0.0083,
     a 0.0052,0.0049,0.0047,0.0067,0.0114,0.0163,0.0124,
     a 0.0062,0.0031,0.0021,0.0021,0.0021,0.0021,0.0021,
     a 0.0021,0.0021,0.0021,0.0021,0.0021,0.0021,0.0021,
     a 0.0021,0.0021,0.0021,0.0031,0.0031,0.0036,0.0049,
     a 0.0061,0.0102,0.0173,0.0336,0.0681,0.1561,0.3157,
     a 0.5213,0.6883,0.7567,0.7837,0.8054,0.8280,0.8446,
     a 0.8582,0.8721,0.8874,0.9052,0.9203,0.9330,0.9450,
     a 0.9546,0.9639,0.9709,0.9761,0.9833,0.9891,0.9949,
     a 0.9999,0.9989,0.9923,0.9797,0.9672,0.9556,0.9446,
     a 0.9371,0.9314,0.9244,0.9177,0.9082,0.8994,0.8935,
     a 0.8854,0.8770,0.8681,0.8630,0.8587,0.8583,0.8577,
     a 0.8590,0.8569,0.8513,0.8271,0.7398,0.5777,0.3807,
     a 0.2127,0.1131,0.0595,0.0315,0.0172,0.0103,0.0066,
     a 0.0043,0.0030,0.0018,0.0010,0.0010,0.0010,0.0010,
     a 0.0010,0.0010,0.0010,0.0010,0.0010,0.0010,0.0010,
     a 0.0010,0.0010,0.0010,0.0010,0.0010,0.0010,0.0010,
     a 0.0010,0.0010,0.0010,0.0010,0.0010,0.0010,0.0010,
     a 0.0010,0.0010,0.0010,0.0010,0.0010,0.0010,0.0010,
     a 0.0010,0.0010,0.0010,0.0010,0.0010,0.0010,0.0010,
     a 0.0010,0.0010,0.0010,0.0010,0.0010,0.0010,0.0010,
     a 0.0010,0.0010,0.0010,0.0010,0.0010,0.0010,0.0010,
     a 0.0010,0.0010,0.0010,0.0010,0.0010,0.0010,0.0010,
     a 0.0010,0.0010,0.0010,0.0010,0.0010,0.0010,0.0010,
     a 0.0010,0.0010,0.0010,0.0010,0.0010,0.0010,0.0010,
     a 0.0010,
     a1160*0./

	wli(1)=0.45
      wls(1)=0.89
      wli(2)=0.45
      wls(2)=0.52
      wli(3)=0.52
      wls(3)=0.59
      wli(4)=0.63
      wls(4)=0.69
      wli(5)=0.77
      wls(5)=0.89
 







      do 1 i=1,1501
      s(i)=sr(iwa,i)
    1 continue
      wlinf=wli(iwa)
      wlsup=wls(iwa)
      return
      end

