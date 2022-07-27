package beast.evolution.substitutionmodel;

import beast.base.core.Description;
import beast.base.evolution.datatype.Aminoacid;
import beast.base.evolution.datatype.DataType;
import beast.base.evolution.substitutionmodel.EmpiricalSubstitutionModel;
/** model data from codonPHYML, which is based on PHYML **/

@Description("LG substitution model for amino acids")
public class OBAMA_LG extends EmpiricalSubstitutionModel {
    @Override
    public double[][] getEmpiricalRates() {
        double[][] rate = new double[20][20];

    
    /* LG model 
     * Si Quang Le  Olivier Gascuel
     * 'An improved general amino-acid replacement matrix' 
     */
    
    rate[1][0] =   0.425093;  rate[2][0] =   0.276818;  rate[2][1] =   0.751878;  rate[3][0] =   0.395144;  
    rate[3][1] =   0.123954;  rate[3][2] =   5.076149;  rate[4][0] =   2.489084;  rate[4][1] =   0.534551;  
    rate[4][2] =   0.528768;  rate[4][3] =   0.062556;  rate[5][0] =   0.969894;  rate[5][1] =   2.807908;  
    rate[5][2] =   1.695752;  rate[5][3] =   0.523386;  rate[5][4] =   0.084808;  rate[6][0] =   1.038545;  
    rate[6][1] =   0.363970;  rate[6][2] =   0.541712;  rate[6][3] =   5.243870;  rate[6][4] =   0.003499;  
    rate[6][5] =   4.128591;  rate[7][0] =   2.066040;  rate[7][1] =   0.390192;  rate[7][2] =   1.437645;  
    rate[7][3] =   0.844926;  rate[7][4] =   0.569265;  rate[7][5] =   0.267959;  rate[7][6] =   0.348847;  
    rate[8][0] =   0.358858;  rate[8][1] =   2.426601;  rate[8][2] =   4.509238;  rate[8][3] =   0.927114;  
    rate[8][4] =   0.640543;  rate[8][5] =   4.813505;  rate[8][6] =   0.423881;  rate[8][7] =   0.311484;  
    rate[9][0] =   0.149830;  rate[9][1] =   0.126991;  rate[9][2] =   0.191503;  rate[9][3] =   0.010690;  
    rate[9][4] =   0.320627;  rate[9][5] =   0.072854;  rate[9][6] =   0.044265;  rate[9][7] =   0.008705;  
    rate[9][8] =   0.108882;  rate[10][0] =   0.395337;  rate[10][1] =   0.301848;  rate[10][2] =   0.068427;  
    rate[10][3] =   0.015076;  rate[10][4] =   0.594007;  rate[10][5] =   0.582457;  rate[10][6] =   0.069673;  
    rate[10][7] =   0.044261;  rate[10][8] =   0.366317;  rate[10][9] =   4.145067;  rate[11][0] =   0.536518;  
    rate[11][1] =   6.326067;  rate[11][2] =   2.145078;  rate[11][3] =   0.282959;  rate[11][4] =   0.013266;  
    rate[11][5] =   3.234294;  rate[11][6] =   1.807177;  rate[11][7] =   0.296636;  rate[11][8] =   0.697264;  
    rate[11][9] =   0.159069;  rate[11][10] =   0.137500;  rate[12][0] =   1.124035;  rate[12][1] =   0.484133;  
    rate[12][2] =   0.371004;  rate[12][3] =   0.025548;  rate[12][4] =   0.893680;  rate[12][5] =   1.672569;  
    rate[12][6] =   0.173735;  rate[12][7] =   0.139538;  rate[12][8] =   0.442472;  rate[12][9] =   4.273607;  
    rate[12][10] =   6.312358;  rate[12][11] =   0.656604;  rate[13][0] =   0.253701;  rate[13][1] =   0.052722;  
    rate[13][2] =   0.089525;  rate[13][3] =   0.017416;  rate[13][4] =   1.105251;  rate[13][5] =   0.035855;  
    rate[13][6] =   0.018811;  rate[13][7] =   0.089586;  rate[13][8] =   0.682139;  rate[13][9] =   1.112727;  
    rate[13][10] =   2.592692;  rate[13][11] =   0.023918;  rate[13][12] =   1.798853;  rate[14][0] =   1.177651;  
    rate[14][1] =   0.332533;  rate[14][2] =   0.161787;  rate[14][3] =   0.394456;  rate[14][4] =   0.075382;  
    rate[14][5] =   0.624294;  rate[14][6] =   0.419409;  rate[14][7] =   0.196961;  rate[14][8] =   0.508851;  
    rate[14][9] =   0.078281;  rate[14][10] =   0.249060;  rate[14][11] =   0.390322;  rate[14][12] =   0.099849;  
    rate[14][13] =   0.094464;  rate[15][0] =   4.727182;  rate[15][1] =   0.858151;  rate[15][2] =   4.008358;  
    rate[15][3] =   1.240275;  rate[15][4] =   2.784478;  rate[15][5] =   1.223828;  rate[15][6] =   0.611973;  
    rate[15][7] =   1.739990;  rate[15][8] =   0.990012;  rate[15][9] =   0.064105;  rate[15][10] =   0.182287;  
    rate[15][11] =   0.748683;  rate[15][12] =   0.346960;  rate[15][13] =   0.361819;  rate[15][14] =   1.338132;  
    rate[16][0] =   2.139501;  rate[16][1] =   0.578987;  rate[16][2] =   2.000679;  rate[16][3] =   0.425860;  
    rate[16][4] =   1.143480;  rate[16][5] =   1.080136;  rate[16][6] =   0.604545;  rate[16][7] =   0.129836;  
    rate[16][8] =   0.584262;  rate[16][9] =   1.033739;  rate[16][10] =   0.302936;  rate[16][11] =   1.136863;  
    rate[16][12] =   2.020366;  rate[16][13] =   0.165001;  rate[16][14] =   0.571468;  rate[16][15] =   6.472279;  
    rate[17][0] =   0.180717;  rate[17][1] =   0.593607;  rate[17][2] =   0.045376;  rate[17][3] =   0.029890;  
    rate[17][4] =   0.670128;  rate[17][5] =   0.236199;  rate[17][6] =   0.077852;  rate[17][7] =   0.268491;  
    rate[17][8] =   0.597054;  rate[17][9] =   0.111660;  rate[17][10] =   0.619632;  rate[17][11] =   0.049906;  
    rate[17][12] =   0.696175;  rate[17][13] =   2.457121;  rate[17][14] =   0.095131;  rate[17][15] =   0.248862;  
    rate[17][16] =   0.140825;  rate[18][0] =   0.218959;  rate[18][1] =   0.314440;  rate[18][2] =   0.612025;  
    rate[18][3] =   0.135107;  rate[18][4] =   1.165532;  rate[18][5] =   0.257336;  rate[18][6] =   0.120037;  
    rate[18][7] =   0.054679;  rate[18][8] =   5.306834;  rate[18][9] =   0.232523;  rate[18][10] =   0.299648;  
    rate[18][11] =   0.131932;  rate[18][12] =   0.481306;  rate[18][13] =   7.803902;  rate[18][14] =   0.089613;  
    rate[18][15] =   0.400547;  rate[18][16] =   0.245841;  rate[18][17] =   3.151815;  rate[19][0] =   2.547870;  
    rate[19][1] =   0.170887;  rate[19][2] =   0.083688;  rate[19][3] =   0.037967;  rate[19][4] =   1.959291;  
    rate[19][5] =   0.210332;  rate[19][6] =   0.245034;  rate[19][7] =   0.076701;  rate[19][8] =   0.119013;  
    rate[19][9] =  10.649107;  rate[19][10] =   1.702745;  rate[19][11] =   0.185202;  rate[19][12] =   1.898718;  
    rate[19][13] =   0.654683;  rate[19][14] =   0.296501;  rate[19][15] =   0.098369;  rate[19][16] =   2.188158;  
    rate[19][17] =   0.189510;  rate[19][18] =   0.249313;  
    
    
       for (int i = 0; i < 20; i++) for (int j = i + 1; j < 20; j++) rate[i][j] = rate[j][i]; return rate;
    }

    @Override
    public double[] getEmpiricalFrequencies() {
        double[] f = new double[20];
    
    f[0] = 0.079066; f[1] = 0.055941; f[2] = 0.041977; f[3] = 0.053052; 
    f[4] = 0.012937; f[5] = 0.040767; f[6] = 0.071586; f[7] = 0.057337; 
    f[8] = 0.022355; f[9] = 0.062157; f[10] = 0.099081; f[11] = 0.064600; 
    f[12] = 0.022951; f[13] = 0.042302; f[14] = 0.044040; f[15] = 0.061197; 
    f[16] = 0.053287; f[17] = 0.012066; f[18] = 0.034155; f[19] = 0.069147; 
    
        return f;
    }

    @Override
    public int[] getEncodingOrder() {
        Aminoacid dataType = new Aminoacid();
        String codeMap = dataType.getCodeMap();
        int[] codeMapNrs = new int[dataType.getStateCount()];
        String encoding = "ARNDCQEGHILKMFPSTWYV";
        for (int i = 0; i < dataType.getStateCount(); i++) {
            codeMapNrs[i] = encoding.indexOf(codeMap.charAt(i));
        }
        return codeMapNrs;
    }

    @Override
    public boolean canHandleDataType(DataType dataType) {
        return dataType instanceof Aminoacid;
    }

/*********************************************************/

} 
