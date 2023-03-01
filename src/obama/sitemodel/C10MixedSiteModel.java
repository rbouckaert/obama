package obama.sitemodel;


import java.util.ArrayList;
import java.util.List;

import beast.base.core.Input.Validate;
import beast.base.core.Log;
import beast.base.evolution.datatype.Aminoacid;
import beast.base.evolution.substitutionmodel.EmpiricalSubstitutionModel;
import beast.base.evolution.substitutionmodel.Frequencies;
import beast.base.evolution.substitutionmodel.GeneralSubstitutionModel;
import beast.base.evolution.substitutionmodel.SubstitutionModel;
import beast.base.inference.parameter.RealParameter;

public class C10MixedSiteModel extends MixedSiteModel {

	double [][] freqs = {{0.408257,0.0349388,0.00698709,0.00978467,0.00616043,0.122161,0.00391518,0.0125784,0.00596702,0.0158339,0.00813132,0.00962854,0.0394156,0.00752797,0.0081783,0.168245,0.0658133,0.0604427,0.00187516,0.00415797},
		{0.102776,0.0149663,0.0155944,0.0419667,0.0180729,0.0138806,0.0158865,0.106608,0.0436344,0.113194,0.04378,0.0213272,0.0223251,0.0440685,0.0418664,0.0529608,0.108174,0.160665,0.00451472,0.0137375},
		{0.0351766,0.00787065,0.000676874,0.00196868,0.0126221,0.00224206,0.00128783,0.351582,0.00188565,0.127818,0.0242632,0.00165915,0.00297716,0.00165596,0.00196786,0.00499981,0.0255378,0.388864,0.00119078,0.00375393},
		{0.0408514,0.00376029,0.233381,0.0901239,0.00251082,0.115833,0.0373197,0.00255236,0.0485017,0.00521646,0.00225718,0.218565,0.0108334,0.0380451,0.0269887,0.0804527,0.030288,0.00444811,0.00108153,0.00698909},
		{0.0185493,0.00704165,0.000977506,0.00248916,0.073333,0.00289529,0.0040104,0.163242,0.00435709,0.444308,0.120282,0.00248957,0.00488276,0.00835394,0.00623624,0.00516424,0.0131807,0.0968581,0.00687598,0.0144734},
		{0.110675,0.00148349,0.163644,0.263846,0.00232568,0.0325228,0.0163804,0.00683349,0.0677158,0.014068,0.00489881,0.0405186,0.0298982,0.0877962,0.035219,0.0562888,0.0426922,0.0181079,0.0010339,0.00405223},
		{0.0522658,0.0143325,0.0297745,0.0388387,0.0624033,0.0228101,0.155164,0.0187406,0.0439469,0.065378,0.0207189,0.0714837,0.0145475,0.073654,0.0668295,0.0549018,0.037014,0.0267512,0.0193757,0.111069},
		{0.0116587,0.0105341,0.00217425,0.00242511,0.365099,0.00347091,0.0366787,0.0187185,0.00266947,0.067649,0.0143535,0.00640111,0.00311599,0.00402037,0.00509901,0.00948485,0.00737139,0.0206341,0.0509565,0.357486},
		{0.0627196,0.00526629,0.0236193,0.0686285,0.00391818,0.0256175,0.0332612,0.0128968,0.227084,0.0305628,0.0124037,0.0428629,0.0140441,0.109811,0.203878,0.0483152,0.0463378,0.0197063,0.00251435,0.00655211},
		{0.114552,0.00985495,0.0416192,0.0364908,0.0046606,0.0503818,0.0165233,0.00929495,0.0423027,0.0139154,0.00822408,0.0750615,0.0379222,0.0339625,0.0324009,0.261065,0.184583,0.0195769,0.0017549,0.00585383},
		};

	public C10MixedSiteModel() {
		mixtureComponentInput.setRule(Validate.FORBIDDEN);
		substModelInput.setRule(Validate.OPTIONAL);
	}

	@Override
    protected List<SubstitutionModel> initialiseMixtureComponents() {
    	List<SubstitutionModel> mixtureComponent = new ArrayList<>();
		Double [] rates = new Double[20*19];
    	if (substModelInput.get() == null) {
    		// assume equal rates
    		Log.warning("Using equal rates for rate matrix");
    		for (int i = 0; i < rates.length; i++) {
    			rates[i] = 1.0;
    		}
    	} else {
    		Log.warning("Using " + substModelInput.get().getClass().getName() + " as rate matrix");
    		double [] r = substModelInput.get().getRateMatrix(null);
    		int [] order = ((EmpiricalSubstitutionModel)substModelInput.get()).getEncodingOrder();
    		int k = 0;
    		for (int i = 0; i < 20; i++) {
        		for (int j = 0; j < 20; j++) {
        			if (i!=j) {
        				rates[k++] = r[order[i]*20+order[j]];
        			}
        		}
    		}
    	}
		int [] order = getEncodingOrder();
		for (int i = 0; i < freqs.length; i++) {
			GeneralSubstitutionModel m = new GeneralSubstitutionModel();
			Double [] f = new Double[20];
			
			for (int j = 0; j < 20; j++) {
				f[j] = this.freqs[i][order[j]];
			}
			
			// normalise to get rid of errors due to frequencies not summing to 1.0 
			// because of precision of freqs values
			double sum = 0;
			for (double d : f) {
				sum +=d;
			}
			for (int j = 0; j < 20; j++) {
				f[j] /= sum;
			}
			
			RealParameter freqs = new RealParameter(f);
			Frequencies frequencies = new Frequencies();
			frequencies.initByName("frequencies", freqs);
			m.initByName("rates", new RealParameter(rates), "frequencies", frequencies);
			mixtureComponent.add(m);
		}
		return mixtureComponent;
	}

	public int[] getEncodingOrder() {
        Aminoacid dataType = new Aminoacid();
        String codeMap = dataType.getCodeMap();
        int[] codeMapNrs = new int[dataType.getStateCount()];
        String encoding = "ACDEFGHIKLMNPQRSTVWY";
        //	String encoding = "ARNDCQEGHILKMFPSTWYV";
        for (int i = 0; i < dataType.getStateCount(); i++) {
            codeMapNrs[i] = encoding.indexOf(codeMap.charAt(i));
        }
        return codeMapNrs;
    }
}
