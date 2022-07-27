package obama.substitutionmodel;

import beast.base.core.Description;
import beast.base.evolution.datatype.Aminoacid;
import beast.base.evolution.datatype.DataType;
import beast.base.evolution.substitutionmodel.EmpiricalSubstitutionModel;
/** model data from codonPHYML, which is based on PHYML **/

@Description("DCMut substitution model for amino acids")
public class OBAMA_DCMut extends EmpiricalSubstitutionModel { 
    @Override
    public double[][] getEmpiricalRates() {
        double[][] rate = new double[20][20];

    /* 
     DCMut : new implementation based on Dayhoff et al.'s raw data and amino acid mutabilities
     C. Kosiol and N. Goldman. (2005),
     ``Different versions of the Dayhoff rate matrix'',
     Mol. Biol. Evol., 22. 193-199.
     */
    rate[1][0] =   26.78280; rate[2][0] =   98.44740; rate[2][1] =   32.70590; rate[3][0] =  119.98050; 
    rate[3][1] =    0.00000; rate[3][2] =  893.15150; rate[4][0] =   36.00160; rate[4][1] =   23.23740; 
    rate[4][2] =    0.00000; rate[4][3] =    0.00000; rate[5][0] =   88.77530; rate[5][1] =  243.99390; 
    rate[5][2] =  102.85090; rate[5][3] =  134.85510; rate[5][4] =    0.00000; rate[6][0] =  196.11670; 
    rate[6][1] =    0.00000; rate[6][2] =  149.34090; rate[6][3] = 1138.86590; rate[6][4] =    0.00000; 
    rate[6][5] =  708.60220; rate[7][0] =  238.61110; rate[7][1] =    8.77910; rate[7][2] =  138.53520; 
    rate[7][3] =  124.09810; rate[7][4] =   10.72780; rate[7][5] =   28.15810; rate[7][6] =   81.19070; 
    rate[8][0] =   22.81160; rate[8][1] =  238.31480; rate[8][2] =  529.00240; rate[8][3] =   86.82410; 
    rate[8][4] =   28.27290; rate[8][5] =  601.16130; rate[8][6] =   43.94690; rate[8][7] =   10.68020; 
    rate[9][0] =   65.34160; rate[9][1] =   63.26290; rate[9][2] =   76.80240; rate[9][3] =   23.92480; 
    rate[9][4] =   43.80740; rate[9][5] =   18.03930; rate[9][6] =   60.95260; rate[9][7] =    0.00000; 
    rate[9][8] =    7.69810; rate[10][0] =   40.64310; rate[10][1] =   15.49240; rate[10][2] =   34.11130; 
    rate[10][3] =    0.00000; rate[10][4] =    0.00000; rate[10][5] =   73.07720; rate[10][6] =   11.28800; 
    rate[10][7] =    7.15140; rate[10][8] =   44.35040; rate[10][9] =  255.66850; rate[11][0] =   25.86350; 
    rate[11][1] =  461.01240; rate[11][2] =  314.83710; rate[11][3] =   71.69130; rate[11][4] =    0.00000; 
    rate[11][5] =  151.90780; rate[11][6] =   83.00780; rate[11][7] =   26.76830; rate[11][8] =   27.04750; 
    rate[11][9] =   46.08570; rate[11][10] =   18.06290; rate[12][0] =   71.78400; rate[12][1] =   89.63210; 
    rate[12][2] =    0.00000; rate[12][3] =    0.00000; rate[12][4] =    0.00000; rate[12][5] =  112.74990; 
    rate[12][6] =   30.48030; rate[12][7] =   17.03720; rate[12][8] =    0.00000; rate[12][9] =  333.27320; 
    rate[12][10] =  523.01150; rate[12][11] =  241.17390; rate[13][0] =   18.36410; rate[13][1] =   13.69060; 
    rate[13][2] =   13.85030; rate[13][3] =    0.00000; rate[13][4] =    0.00000; rate[13][5] =    0.00000; 
    rate[13][6] =    0.00000; rate[13][7] =   15.34780; rate[13][8] =   47.59270; rate[13][9] =  195.19510; 
    rate[13][10] =  156.51600; rate[13][11] =    0.00000; rate[13][12] =   92.18600; rate[14][0] =  248.59200; 
    rate[14][1] =  102.83130; rate[14][2] =   41.92440; rate[14][3] =   13.39400; rate[14][4] =   18.75500; 
    rate[14][5] =  152.61880; rate[14][6] =   50.70030; rate[14][7] =   34.71530; rate[14][8] =   93.37090; 
    rate[14][9] =   11.91520; rate[14][10] =   31.62580; rate[14][11] =   33.54190; rate[14][12] =   17.02050; 
    rate[14][13] =   11.05060; rate[15][0] =  405.18700; rate[15][1] =  153.15900; rate[15][2] =  488.58920; 
    rate[15][3] =   95.60970; rate[15][4] =  159.83560; rate[15][5] =   56.18280; rate[15][6] =   79.39990; 
    rate[15][7] =  232.22430; rate[15][8] =   35.36430; rate[15][9] =   24.79550; rate[15][10] =   17.14320; 
    rate[15][11] =   95.45570; rate[15][12] =   61.99510; rate[15][13] =   45.99010; rate[15][14] =  242.72020; 
    rate[16][0] =  368.03650; rate[16][1] =   26.57450; rate[16][2] =  227.16970; rate[16][3] =   66.09300; 
    rate[16][4] =   16.23660; rate[16][5] =   52.56510; rate[16][6] =   34.01560; rate[16][7] =   30.66620; 
    rate[16][8] =   22.63330; rate[16][9] =  190.07390; rate[16][10] =   33.10900; rate[16][11] =  135.05990; 
    rate[16][12] =  103.15340; rate[16][13] =   13.66550; rate[16][14] =   78.28570; rate[16][15] =  543.66740; 
    rate[17][0] =    0.00000; rate[17][1] =  200.13750; rate[17][2] =   22.49680; rate[17][3] =    0.00000; 
    rate[17][4] =    0.00000; rate[17][5] =    0.00000; rate[17][6] =    0.00000; rate[17][7] =    0.00000; 
    rate[17][8] =   27.05640; rate[17][9] =    0.00000; rate[17][10] =   46.17760; rate[17][11] =    0.00000; 
    rate[17][12] =    0.00000; rate[17][13] =   76.23540; rate[17][14] =    0.00000; rate[17][15] =   74.08190; 
    rate[17][16] =    0.00000; rate[18][0] =   24.41390; rate[18][1] =    7.80120; rate[18][2] =   94.69400; 
    rate[18][3] =    0.00000; rate[18][4] =   95.31640; rate[18][5] =    0.00000; rate[18][6] =   21.47170; 
    rate[18][7] =    0.00000; rate[18][8] =  126.54000; rate[18][9] =   37.48340; rate[18][10] =   28.65720; 
    rate[18][11] =   13.21420; rate[18][12] =    0.00000; rate[18][13] =  695.26290; rate[18][14] =    0.00000; 
    rate[18][15] =   33.62890; rate[18][16] =   41.78390; rate[18][17] =   60.80700; rate[19][0] =  205.95640; 
    rate[19][1] =   24.03680; rate[19][2] =   15.80670; rate[19][3] =   17.83160; rate[19][4] =   48.46780; 
    rate[19][5] =   34.69830; rate[19][6] =   36.72500; rate[19][7] =   53.81650; rate[19][8] =   43.87150; 
    rate[19][9] =  881.00380; rate[19][10] =  174.51560; rate[19][11] =   10.38500; rate[19][12] =  256.59550; 
    rate[19][13] =   12.36060; rate[19][14] =   48.50260; rate[19][15] =   30.38360; rate[19][16] =  156.19970; 
    rate[19][17] =    0.00000; rate[19][18] =   27.93790; 
    
       for (int i = 0; i < 20; i++) for (int j = i + 1; j < 20; j++) rate[i][j] = rate[j][i]; return rate;
    }

    @Override
    public double[] getEmpiricalFrequencies() {
        double[] f = new double[20];
    
    f[ 0] = 0.087127; f[ 1] = 0.040904; f[ 2] = 0.040432; f[ 3] = 0.046872; 
    f[ 4] = 0.033474; f[ 5] = 0.038255; f[ 6] = 0.049530; f[ 7] = 0.088612;
    f[ 8] = 0.033619; f[ 9] = 0.036886; f[10] = 0.085357; f[11] = 0.080481;
    f[12] = 0.014753; f[13] = 0.039772; f[14] = 0.050680; f[15] = 0.069577;
    f[16] = 0.058542; f[17] = 0.010494; f[18] = 0.029916; f[19] = 0.064717; // was 0.064718, but does not add to 1
    
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

} 
