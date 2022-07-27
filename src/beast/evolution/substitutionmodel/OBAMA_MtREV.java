package beast.evolution.substitutionmodel;

import beast.base.core.Description;
import beast.base.evolution.datatype.Aminoacid;
import beast.base.evolution.datatype.DataType;
import beast.base.evolution.substitutionmodel.EmpiricalSubstitutionModel;
/** model data from codonPHYML, which is based on PHYML **/

@Description("MtREV substitution model for amino acids")
public class OBAMA_MtREV extends EmpiricalSubstitutionModel {
    @Override
    public double[][] getEmpiricalRates() {
        double[][] rate = new double[20][20];

    /* J. Adachi and M. Hasegawa, ``Model of amino acid substitution in proteins
     encoded by mitochondrial DNA'' J. Mol. Evol. 42, 459 (1996) */
    
    
    
    rate[1][0] =   23.18; rate[2][0] =   26.95; rate[2][1] =   13.24; rate[3][0] =   17.67;
    rate[3][1] =    1.90; rate[3][2] =  794.38; rate[4][0] =   59.93; rate[4][1] =  103.33;
    rate[4][2] =   58.94; rate[4][3] =    1.90; rate[5][0] =    1.90; rate[5][1] =  220.99;
    rate[5][2] =  173.56; rate[5][3] =   55.28; rate[5][4] =   75.24; rate[6][0] =    9.77;
    rate[6][1] =    1.90; rate[6][2] =   63.05; rate[6][3] =  583.55; rate[6][4] =    1.90;
    rate[6][5] =  313.56; rate[7][0] =  120.71; rate[7][1] =   23.03; rate[7][2] =   53.30;
    rate[7][3] =   56.77; rate[7][4] =   30.71; rate[7][5] =    6.75; rate[7][6] =   28.28;
    rate[8][0] =   13.90; rate[8][1] =  165.23; rate[8][2] =  496.13; rate[8][3] =  113.99;
    rate[8][4] =  141.49; rate[8][5] =  582.40; rate[8][6] =   49.12; rate[8][7] =    1.90;
    rate[9][0] =   96.49; rate[9][1] =    1.90; rate[9][2] =   27.10; rate[9][3] =    4.34;
    rate[9][4] =   62.73; rate[9][5] =    8.34; rate[9][6] =    3.31; rate[9][7] =    5.98;
    rate[9][8] =   12.26; rate[10][0] =   25.46; rate[10][1] =   15.58; rate[10][2] =   15.16;
    rate[10][3] =    1.90; rate[10][4] =   25.65; rate[10][5] =   39.70; rate[10][6] =    1.90;
    rate[10][7] =    2.41; rate[10][8] =   11.49; rate[10][9] =  329.09; rate[11][0] =    8.36;
    rate[11][1] =  141.40; rate[11][2] =  608.70; rate[11][3] =    2.31; rate[11][4] =    1.90;
    rate[11][5] =  465.58; rate[11][6] =  313.86; rate[11][7] =   22.73; rate[11][8] =  127.67;
    rate[11][9] =   19.57; rate[11][10] =   14.88; rate[12][0] =  141.88; rate[12][1] =    1.90;
    rate[12][2] =   65.41; rate[12][3] =    1.90; rate[12][4] =    6.18; rate[12][5] =   47.37;
    rate[12][6] =    1.90; rate[12][7] =    1.90; rate[12][8] =   11.97; rate[12][9] =  517.98;
    rate[12][10] =  537.53; rate[12][11] =   91.37; rate[13][0] =    6.37; rate[13][1] =    4.69;
    rate[13][2] =   15.20; rate[13][3] =    4.98; rate[13][4] =   70.80; rate[13][5] =   19.11;
    rate[13][6] =    2.67; rate[13][7] =    1.90; rate[13][8] =   48.16; rate[13][9] =   84.67;
    rate[13][10] =  216.06; rate[13][11] =    6.44; rate[13][12] =   90.82; rate[14][0] =   54.31;
    rate[14][1] =   23.64; rate[14][2] =   73.31; rate[14][3] =   13.43; rate[14][4] =   31.26;
    rate[14][5] =  137.29; rate[14][6] =   12.83; rate[14][7] =    1.90; rate[14][8] =   60.97;
    rate[14][9] =   20.63; rate[14][10] =   40.10; rate[14][11] =   50.10; rate[14][12] =   18.84;
    rate[14][13] =   17.31; rate[15][0] =  387.86; rate[15][1] =    6.04; rate[15][2] =  494.39;
    rate[15][3] =   69.02; rate[15][4] =  277.05; rate[15][5] =   54.11; rate[15][6] =   54.71;
    rate[15][7] =  125.93; rate[15][8] =   77.46; rate[15][9] =   47.70; rate[15][10] =   73.61;
    rate[15][11] =  105.79; rate[15][12] =  111.16; rate[15][13] =   64.29; rate[15][14] =  169.90;
    rate[16][0] =  480.72; rate[16][1] =    2.08; rate[16][2] =  238.46; rate[16][3] =   28.01;
    rate[16][4] =  179.97; rate[16][5] =   94.93; rate[16][6] =   14.82; rate[16][7] =   11.17;
    rate[16][8] =   44.78; rate[16][9] =  368.43; rate[16][10] =  126.40; rate[16][11] =  136.33;
    rate[16][12] =  528.17; rate[16][13] =   33.85; rate[16][14] =  128.22; rate[16][15] =  597.21;
    rate[17][0] =    1.90; rate[17][1] =   21.95; rate[17][2] =   10.68; rate[17][3] =   19.86;
    rate[17][4] =   33.60; rate[17][5] =    1.90; rate[17][6] =    1.90; rate[17][7] =   10.92;
    rate[17][8] =    7.08; rate[17][9] =    1.90; rate[17][10] =   32.44; rate[17][11] =   24.00;
    rate[17][12] =   21.71; rate[17][13] =    7.84; rate[17][14] =    4.21; rate[17][15] =   38.58;
    rate[17][16] =    9.99; rate[18][0] =    6.48; rate[18][1] =    1.90; rate[18][2] =  191.36;
    rate[18][3] =   21.21; rate[18][4] =  254.77; rate[18][5] =   38.82; rate[18][6] =   13.12;
    rate[18][7] =    3.21; rate[18][8] =  670.14; rate[18][9] =   25.01; rate[18][10] =   44.15;
    rate[18][11] =   51.17; rate[18][12] =   39.96; rate[18][13] =  465.58; rate[18][14] =   16.21;
    rate[18][15] =   64.92; rate[18][16] =   38.73; rate[18][17] =   26.25; rate[19][0] =  195.06;
    rate[19][1] =    7.64; rate[19][2] =    1.90; rate[19][3] =    1.90; rate[19][4] =    1.90;
    rate[19][5] =   19.00; rate[19][6] =   21.14; rate[19][7] =    2.53; rate[19][8] =    1.90;
    rate[19][9] = 1222.94; rate[19][10] =   91.67; rate[19][11] =    1.90; rate[19][12] =  387.54;
    rate[19][13] =    6.35; rate[19][14] =    8.23; rate[19][15] =    1.90; rate[19][16] =  204.54;
    rate[19][17] =    5.37; rate[19][18] =    1.90;
    
       for (int i = 0; i < 20; i++) for (int j = i + 1; j < 20; j++) rate[i][j] = rate[j][i]; return rate;
    }

    @Override
    public double[] getEmpiricalFrequencies() {
        double[] f = new double[20];
    
    f[ 0] = 0.072000; f[ 1] = 0.019000; f[ 2] = 0.039000; f[ 3] = 0.019000;
    f[ 4] = 0.006000; f[ 5] = 0.025000; f[ 6] = 0.024000; f[ 7] = 0.056000;
    f[ 8] = 0.028000; f[ 9] = 0.088000; f[10] = 0.169000; f[11] = 0.023000;
    f[12] = 0.054000; f[13] = 0.061000; f[14] = 0.054000; f[15] = 0.072000;
    f[16] = 0.086000; f[17] = 0.029000; f[18] = 0.033000; f[19] = 0.043000;
    
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
