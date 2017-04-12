package beast.evolution.substitutionmodel;

import beast.evolution.datatype.Aminoacid;
import beast.evolution.datatype.DataType;
/** model data from codonPHYML, which is based on PHYML **/

public class BAMA_WAG extends EmpiricalSubstitutionModel {
    @Override
    double[][] getEmpiricalRates() {
        double[][] rate = new double[20][20];

    
    /* WAG's model data
     * Simon Whelan and Nick Goldman
     * 'A general empirical model of protein evolution derived from multiple
     *  protein families using a maximum-likelihood approach' 
     * MBE (2001) 18:691-699
     */
    
    
    rate[1][0] =  55.15710; rate[2][0] =  50.98480; rate[2][1] =  63.53460; 
    rate[3][0] =  73.89980; rate[3][1] =  14.73040; rate[3][2] = 542.94200; 
    rate[4][0] = 102.70400; rate[4][1] =  52.81910; rate[4][2] =  26.52560; 
    rate[4][3] =   3.02949; rate[5][0] =  90.85980; rate[5][1] = 303.55000; 
    rate[5][2] = 154.36400; rate[5][3] =  61.67830; rate[5][4] =   9.88179; 
    rate[6][0] = 158.28500; rate[6][1] =  43.91570; rate[6][2] =  94.71980; 
    rate[6][3] = 617.41600; rate[6][4] =   2.13520; rate[6][5] = 546.94700; 
    rate[7][0] = 141.67200; rate[7][1] =  58.46650; rate[7][2] = 112.55600; 
    rate[7][3] =  86.55840; rate[7][4] =  30.66740; rate[7][5] =  33.00520; 
    rate[7][6] =  56.77170; rate[8][0] =  31.69540; rate[8][1] = 213.71500; 
    rate[8][2] = 395.62900; rate[8][3] =  93.06760; rate[8][4] =  24.89720; 
    rate[8][5] = 429.41100; rate[8][6] =  57.00250; rate[8][7] =  24.94100; 
    rate[9][0] =  19.33350; rate[9][1] =  18.69790; rate[9][2] =  55.42360; 
    rate[9][3] =   3.94370; rate[9][4] =  17.01350; rate[9][5] =  11.39170; 
    rate[9][6] =  12.73950; rate[9][7] =   3.04501; rate[9][8] =  13.81900; 
    rate[10][0] =  39.79150; rate[10][1] =  49.76710; rate[10][2] =  13.15280; 
    rate[10][3] =   8.48047; rate[10][4] =  38.42870; rate[10][5] =  86.94890; 
    rate[10][6] =  15.42630; rate[10][7] =   6.13037; rate[10][8] =  49.94620; 
    rate[10][9] = 317.09700; rate[11][0] =  90.62650; rate[11][1] = 535.14200; 
    rate[11][2] = 301.20100; rate[11][3] =  47.98550; rate[11][4] =   7.40339; 
    rate[11][5] = 389.49000; rate[11][6] = 258.44300; rate[11][7] =  37.35580; 
    rate[11][8] =  89.04320; rate[11][9] =  32.38320; rate[11][10] =  25.75550; 
    rate[12][0] =  89.34960; rate[12][1] =  68.31620; rate[12][2] =  19.82210; 
    rate[12][3] =  10.37540; rate[12][4] =  39.04820; rate[12][5] = 154.52600; 
    rate[12][6] =  31.51240; rate[12][7] =  17.41000; rate[12][8] =  40.41410; 
    rate[12][9] = 425.74600; rate[12][10] = 485.40200; rate[12][11] =  93.42760; 
    rate[13][0] =  21.04940; rate[13][1] =  10.27110; rate[13][2] =   9.61621; 
    rate[13][3] =   4.67304; rate[13][4] =  39.80200; rate[13][5] =   9.99208; 
    rate[13][6] =   8.11339; rate[13][7] =   4.99310; rate[13][8] =  67.93710; 
    rate[13][9] = 105.94700; rate[13][10] = 211.51700; rate[13][11] =   8.88360; 
    rate[13][12] = 119.06300; rate[14][0] = 143.85500; rate[14][1] =  67.94890; 
    rate[14][2] =  19.50810; rate[14][3] =  42.39840; rate[14][4] =  10.94040; 
    rate[14][5] =  93.33720; rate[14][6] =  68.23550; rate[14][7] =  24.35700; 
    rate[14][8] =  69.61980; rate[14][9] =   9.99288; rate[14][10] =  41.58440; 
    rate[14][11] =  55.68960; rate[14][12] =  17.13290; rate[14][13] =  16.14440; 
    rate[15][0] = 337.07900; rate[15][1] = 122.41900; rate[15][2] = 397.42300; 
    rate[15][3] = 107.17600; rate[15][4] = 140.76600; rate[15][5] = 102.88700; 
    rate[15][6] =  70.49390; rate[15][7] = 134.18200; rate[15][8] =  74.01690; 
    rate[15][9] =  31.94400; rate[15][10] =  34.47390; rate[15][11] =  96.71300; 
    rate[15][12] =  49.39050; rate[15][13] =  54.59310; rate[15][14] = 161.32800; 
    rate[16][0] = 212.11100; rate[16][1] =  55.44130; rate[16][2] = 203.00600; 
    rate[16][3] =  37.48660; rate[16][4] =  51.29840; rate[16][5] =  85.79280; 
    rate[16][6] =  82.27650; rate[16][7] =  22.58330; rate[16][8] =  47.33070; 
    rate[16][9] = 145.81600; rate[16][10] =  32.66220; rate[16][11] = 138.69800; 
    rate[16][12] = 151.61200; rate[16][13] =  17.19030; rate[16][14] =  79.53840; 
    rate[16][15] = 437.80200; rate[17][0] =  11.31330; rate[17][1] = 116.39200; 
    rate[17][2] =   7.19167; rate[17][3] =  12.97670; rate[17][4] =  71.70700; 
    rate[17][5] =  21.57370; rate[17][6] =  15.65570; rate[17][7] =  33.69830; 
    rate[17][8] =  26.25690; rate[17][9] =  21.24830; rate[17][10] =  66.53090; 
    rate[17][11] =  13.75050; rate[17][12] =  51.57060; rate[17][13] = 152.96400; 
    rate[17][14] =  13.94050; rate[17][15] =  52.37420; rate[17][16] =  11.08640; 
    rate[18][0] =  24.07350; rate[18][1] =  38.15330; rate[18][2] = 108.60000; 
    rate[18][3] =  32.57110; rate[18][4] =  54.38330; rate[18][5] =  22.77100; 
    rate[18][6] =  19.63030; rate[18][7] =  10.36040; rate[18][8] = 387.34400; 
    rate[18][9] =  42.01700; rate[18][10] =  39.86180; rate[18][11] =  13.32640; 
    rate[18][12] =  42.84370; rate[18][13] = 645.42800; rate[18][14] =  21.60460; 
    rate[18][15] =  78.69930; rate[18][16] =  29.11480; rate[18][17] = 248.53900; 
    rate[19][0] = 200.60100; rate[19][1] =  25.18490; rate[19][2] =  19.62460; 
    rate[19][3] =  15.23350; rate[19][4] = 100.21400; rate[19][5] =  30.12810; 
    rate[19][6] =  58.87310; rate[19][7] =  18.72470; rate[19][8] =  11.83580; 
    rate[19][9] = 782.13000; rate[19][10] = 180.03400; rate[19][11] =  30.54340; 
    rate[19][12] = 205.84500; rate[19][13] =  64.98920; rate[19][14] =  31.48870; 
    rate[19][15] =  23.27390; rate[19][16] = 138.82300; rate[19][17] =  36.53690; 
    rate[19][18] =  31.47300; 
    
       for (int i = 0; i < 20; i++) for (int j = i + 1; j < 20; j++) rate[i][j] = rate[j][i]; return rate;
    }

    @Override
    public double[] getEmpiricalFrequencies() {
        double[] f = new double[20];
    
    f[0] = 0.0866279; f[1] =  0.043972; f[2] =  0.0390894; f[3] =  0.0570451;
    f[4] =  0.0193078; f[5] =  0.0367281; f[6] =  0.0580589; f[7] =  0.0832518;
    f[8] =  0.0244313; f[9] =  0.048466; f[10] =  0.086209; f[11] = 0.0620286;
    f[12] = 0.0195027; f[13] =  0.0384319; f[14] =  0.0457631; f[15] = 0.0695179;
    f[16] =  0.0610127; f[17] =  0.0143859; f[18] =  0.0352742; f[19] =  0.0708956;
    
    
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
