package obama.sitemodel;

public class C40MixedSiteModel extends C10MixedSiteModel {

	public C40MixedSiteModel() {
		super();
		freqs = new double[][]{{0.066026,0.00565867,0.105447,0.0440361,0.00131048,0.0711239,0.0168195,0.00390887,0.036669,0.0055316,0.00374124,0.159982,0.0176359,0.0273928,0.0231862,0.249769,0.150708,0.0065529,0.000672321,0.00382902},
			{0.0232377,0.00379875,0.353209,0.0739378,0.00240321,0.0576668,0.0315867,0.00310928,0.0259363,0.00387116,0.00173556,0.275965,0.00631169,0.0197339,0.0122683,0.0657068,0.0270484,0.00475317,0.000760289,0.00696025},
			{0.0166487,0.00366657,0.000565145,0.00133563,0.00827757,0.000889475,0.000823185,0.412937,0.00119041,0.0884689,0.0186055,0.00126222,0.001403,0.00106698,0.00125948,0.00213394,0.0162167,0.420686,0.000608205,0.00195532},
			{0.239474,0.0283812,0.00447417,0.010553,0.00559911,0.013511,0.00389298,0.0765957,0.0071093,0.0358495,0.0199496,0.0120537,0.0114266,0.00865589,0.00729013,0.0847799,0.179728,0.245468,0.0009838,0.00422407},
			{0.119461,0.0150527,0.0134273,0.0192173,0.0550467,0.0337676,0.0214746,0.0579002,0.0147261,0.144631,0.0561243,0.0294552,0.0631355,0.0301538,0.0233256,0.0925267,0.083123,0.0811758,0.0131636,0.0331118},
			{0.0567044,0.00089248,0.29555,0.379515,0.00129723,0.023047,0.0118361,0.0031182,0.0314206,0.00601375,0.00285841,0.0364734,0.0124746,0.0609517,0.0117359,0.0300335,0.0227051,0.00946396,0.000773876,0.00313438},
			{0.0179027,0.016076,0.000887041,0.00231821,0.334486,0.00398298,0.0127293,0.0404651,0.00279947,0.167614,0.0424172,0.00356977,0.00201151,0.00453955,0.00409671,0.00758416,0.00682273,0.0326045,0.0518381,0.245254},
			{0.271217,0.200383,0.0021017,0.002323,0.020299,0.0502501,0.0053728,0.0150685,0.00206463,0.0330003,0.0154811,0.0141045,0.0045351,0.00482641,0.00564808,0.17642,0.0839578,0.0741934,0.00462652,0.0141271},
			{0.0894737,0.00455383,0.0272183,0.127508,0.00565902,0.0115686,0.0215746,0.0469424,0.138205,0.0512035,0.0147657,0.0190192,0.00955465,0.116809,0.104003,0.0383954,0.0836653,0.0819556,0.00170794,0.00621813},
			{0.0495441,0.0182506,0.0143641,0.0215379,0.141805,0.01402,0.110854,0.0247066,0.0258142,0.0700288,0.0188272,0.0315864,0.0112101,0.0316504,0.0375346,0.0456094,0.0361428,0.0369178,0.0371985,0.222397},
			{0.170431,0.000974733,0.109856,0.253646,0.00133213,0.0249846,0.010139,0.00587494,0.0903324,0.0116526,0.00365127,0.0271109,0.0293614,0.09173,0.0415784,0.0561766,0.0479046,0.02033,0.000669682,0.00226373},
			{0.0162725,0.00506141,0.00101821,0.00251413,0.0376246,0.00219354,0.00299143,0.132817,0.00401204,0.490444,0.192993,0.00218762,0.00343332,0.0104414,0.00548261,0.00401221,0.0127074,0.064772,0.00321076,0.00581006},
			{0.0823766,0.00656943,0.0311745,0.0675531,0.00647179,0.0178962,0.0251144,0.0291162,0.0982302,0.0287904,0.0168023,0.059839,0.0114045,0.0686451,0.0734226,0.1303,0.182037,0.0540271,0.00227246,0.00795733},
			{0.359497,0.0251417,0.00314844,0.00649627,0.00920205,0.119468,0.00229704,0.0458767,0.00501688,0.0468054,0.0215569,0.00334215,0.0443916,0.00490143,0.00724072,0.0465271,0.0477755,0.194216,0.00245402,0.00464504},
			{0.201558,0.00323653,0.095415,0.153491,0.00256433,0.0667292,0.0155219,0.00677408,0.0547323,0.0165114,0.0060163,0.0425386,0.00919706,0.0772011,0.0430162,0.118598,0.0625473,0.0202798,0.000956551,0.003115},
			{0.104273,0.00224501,0.242407,0.177482,0.00125703,0.169782,0.0132649,0.00189295,0.0220652,0.00425426,0.00164412,0.0621646,0.0317042,0.0356499,0.0147062,0.0778636,0.0288516,0.00602502,0.00069309,0.00177419},
			{0.0781183,0.0194449,0.00415417,0.0116634,0.0262794,0.0111524,0.00635894,0.135453,0.00937298,0.245757,0.108778,0.015927,0.0055294,0.0240152,0.0111498,0.0408519,0.0860514,0.148276,0.00315476,0.00851085},
			{0.0856592,0.0136073,0.0135062,0.00786026,0.0047153,0.0245401,0.00553791,0.0100592,0.0127319,0.0103344,0.00806758,0.0441923,0.0175274,0.00925906,0.0101233,0.340648,0.357329,0.019367,0.00142431,0.00350998},
			{0.0674595,0.00216342,0.0662588,0.0865501,0.00182127,0.0368557,0.0381149,0.00332388,0.189974,0.009384,0.00394874,0.116311,0.0151208,0.093936,0.116173,0.0842204,0.0565954,0.00645142,0.00071873,0.00461894},
			{0.0572262,0.00153015,0.179393,0.199226,0.00137018,0.0316472,0.0291392,0.00458046,0.101562,0.010074,0.00402046,0.108388,0.00636741,0.0903669,0.0494724,0.0621143,0.0496102,0.00859413,0.000666929,0.00464976},
			{0.00360202,0.00454848,0.00208716,0.00178577,0.0855715,0.00563916,0.00649688,0.00292929,0.00104198,0.0232635,0.00445923,0.00134555,0.0024992,0.00327181,0.0102713,0.00306718,0.00259003,0.00586684,0.761782,0.067881},
			{0.203202,0.00981316,0.0135012,0.00838182,0.00196196,0.618489,0.00277479,0.00118285,0.00445989,0.00398268,0.00206318,0.0143744,0.00858704,0.00445146,0.00838957,0.073992,0.0108922,0.00607691,0.00186061,0.00156387},
			{0.00508988,0.00617765,0.00161262,0.00120404,0.356313,0.00163342,0.0393461,0.00590888,0.00137137,0.0249344,0.00497952,0.0057093,0.00141364,0.00246931,0.00287408,0.00595277,0.00365368,0.00819341,0.0357987,0.485365},
			{0.0403336,0.00815495,0.00982186,0.0375407,0.0119141,0.00479344,0.0176736,0.189342,0.0607377,0.105186,0.03056,0.0216052,0.00775506,0.0383639,0.0540186,0.025711,0.100991,0.221091,0.002878,0.0115277},
			{0.0790086,0.00249479,0.0546012,0.199788,0.00237734,0.0192656,0.02707,0.00756675,0.155311,0.0254542,0.00980244,0.0309384,0.00566407,0.184338,0.106544,0.0332371,0.0359575,0.0145306,0.00116828,0.00488208},
			{0.0722241,0.00647553,0.119488,0.134589,0.0348233,0.0287815,0.0699011,0.0173589,0.0490342,0.051987,0.0154411,0.067893,0.0145597,0.070897,0.0489728,0.058958,0.0425973,0.0317884,0.00879138,0.0554387},
			{0.108584,0.00125026,0.152967,0.166485,0.00145535,0.0336098,0.0134902,0.00388218,0.0576227,0.00898614,0.0024339,0.0441956,0.19901,0.0405398,0.020645,0.084675,0.0454715,0.0113416,0.000590283,0.00276502},
			{0.0309526,0.0065594,0.00823521,0.0291974,0.00368916,0.0154206,0.0310385,0.00982516,0.306263,0.02379,0.00970717,0.0301337,0.00950291,0.0832608,0.319589,0.0295285,0.0303052,0.0133037,0.00281253,0.00688506},
			{0.00989537,0.00282767,0.000374823,0.00091821,0.0298607,0.000699707,0.00104195,0.311504,0.00139605,0.375039,0.0474451,0.000730793,0.00252963,0.0017337,0.00196045,0.0014628,0.0075739,0.1973,0.00167998,0.00402599},
			{0.116321,0.00347923,0.0731918,0.138088,0.00941177,0.0193193,0.0160241,0.0712243,0.035512,0.0771474,0.0242841,0.0250164,0.0508927,0.0586677,0.0273321,0.047556,0.0726552,0.123571,0.00268927,0.0076166},
			{0.128522,0.0172929,0.040275,0.0250692,0.0199118,0.112703,0.0606981,0.010935,0.028875,0.0258416,0.0167593,0.117984,0.0180675,0.0439706,0.0373073,0.174149,0.0648968,0.0182067,0.0063575,0.0321772},
			{0.0372287,0.014494,0.00237032,0.00485851,0.146377,0.00464339,0.0186795,0.182046,0.00581985,0.17801,0.0371334,0.00533773,0.00485386,0.00790971,0.0094528,0.0103571,0.0284162,0.184992,0.0211294,0.0958905},
			{0.0535644,0.00962562,0.0113537,0.0391699,0.0264214,0.0120279,0.0384888,0.0522748,0.0996038,0.189239,0.0712219,0.0239173,0.00837206,0.0928585,0.11598,0.0299114,0.0389485,0.0500948,0.0104232,0.026503},
			{0.133242,0.0246514,0.00127392,0.00404615,0.0156995,0.00891392,0.00158647,0.197128,0.00237132,0.125129,0.0286947,0.0022705,0.0118846,0.00308435,0.00331477,0.0171462,0.0563298,0.356621,0.00173418,0.00487784},
			{0.149866,0.00354374,0.0280355,0.0435381,0.00475339,0.0311113,0.0140626,0.0101953,0.0393125,0.0251434,0.00515483,0.0176453,0.39238,0.0348151,0.0326607,0.0874497,0.0473307,0.0271597,0.00152152,0.00432083},
			{0.421437,0.018761,0.00733051,0.00868378,0.00271561,0.0902333,0.0030262,0.00393628,0.00515087,0.00471933,0.00383066,0.012159,0.020894,0.00727486,0.0061426,0.290119,0.0651922,0.0252211,0.000810824,0.00236228},
			{0.177071,0.00783489,0.02265,0.0509767,0.00405142,0.089739,0.0220667,0.00595198,0.125769,0.020537,0.00929825,0.0311657,0.0264088,0.0752471,0.133278,0.116959,0.0565567,0.0165087,0.00299471,0.00493467},
			{0.0293984,0.00317291,0.109971,0.046427,0.00150396,0.422242,0.0272495,0.000799733,0.0622314,0.00376343,0.00166571,0.148362,0.00564818,0.0388688,0.0370902,0.0472252,0.0086569,0.00203639,0.000917602,0.00276931},
			{0.0265779,0.0101369,0.0280314,0.0269057,0.0276961,0.0173377,0.281513,0.0064647,0.0474749,0.026821,0.00723753,0.13186,0.0083015,0.0989711,0.0791105,0.0426277,0.0259043,0.0100147,0.00785289,0.0891598},
			{0.00960965,0.00773017,0.000633186,0.00104719,0.263017,0.00202274,0.00390014,0.0733098,0.00149315,0.445169,0.0732575,0.00131044,0.00427681,0.00338994,0.00271362,0.00361174,0.00579284,0.0425173,0.0181276,0.0370698},
			};
	}

}
