package obama.substitutionmodel;

import beast.base.core.Description;
import beast.base.evolution.datatype.Aminoacid;
import beast.base.evolution.datatype.DataType;
import beast.base.evolution.substitutionmodel.EmpiricalSubstitutionModel;
/** model data from codonPHYML, which is based on PHYML **/

@Description("Dayhoff substitution model for amino acids")
public class OBAMA_Dayhoff extends EmpiricalSubstitutionModel { 
    @Override
    public double[][] getEmpiricalRates() {
        double[][] rate = new double[20][20];

    /* Dayhoff's model data
     * Dayhoff, M.O., Schwartz, R.M., Orcutt, B.C. (1978)
     * "A model of evolutionary change in proteins."
     * Dayhoff, M.O.(ed.) Atlas of Protein Sequence Structur., Vol5, Suppl3.
     * National Biomedical Research Foundation, Washington DC, pp.345-352.
     */
    rate[1][0] =   27.00; rate[2][0] =   98.00; rate[2][1] =   32.00; rate[3][0] =  120.00;
    rate[3][1] =    0.00; rate[3][2] =  905.00; rate[4][0] =   36.00; rate[4][1] =   23.00;
    rate[4][2] =    0.00; rate[4][3] =    0.00; rate[5][0] =   89.00; rate[5][1] =  246.00;
    rate[5][2] =  103.00; rate[5][3] =  134.00; rate[5][4] =    0.00; rate[6][0] =  198.00;
    rate[6][1] =    1.00; rate[6][2] =  148.00; rate[6][3] = 1153.00; rate[6][4] =    0.00;
    rate[6][5] =  716.00; rate[7][0] =  240.00; rate[7][1] =    9.00; rate[7][2] =  139.00;
    rate[7][3] =  125.00; rate[7][4] =   11.00; rate[7][5] =   28.00; rate[7][6] =   81.00;
    rate[8][0] =   23.00; rate[8][1] =  240.00; rate[8][2] =  535.00; rate[8][3] =   86.00;
    rate[8][4] =   28.00; rate[8][5] =  606.00; rate[8][6] =   43.00; rate[8][7] =   10.00;
    rate[9][0] =   65.00; rate[9][1] =   64.00; rate[9][2] =   77.00; rate[9][3] =   24.00;
    rate[9][4] =   44.00; rate[9][5] =   18.00; rate[9][6] =   61.00; rate[9][7] =    0.00;
    rate[9][8] =    7.00; rate[10][0] =   41.00; rate[10][1] =   15.00; rate[10][2] =   34.00;
    rate[10][3] =    0.00; rate[10][4] =    0.00; rate[10][5] =   73.00; rate[10][6] =   11.00;
    rate[10][7] =    7.00; rate[10][8] =   44.00; rate[10][9] =  257.00; rate[11][0] =   26.00;
    rate[11][1] =  464.00; rate[11][2] =  318.00; rate[11][3] =   71.00; rate[11][4] =    0.00;
    rate[11][5] =  153.00; rate[11][6] =   83.00; rate[11][7] =   27.00; rate[11][8] =   26.00;
    rate[11][9] =   46.00; rate[11][10] =   18.00; rate[12][0] =   72.00; rate[12][1] =   90.00;
    rate[12][2] =    1.00; rate[12][3] =    0.00; rate[12][4] =    0.00; rate[12][5] =  114.00;
    rate[12][6] =   30.00; rate[12][7] =   17.00; rate[12][8] =    0.00; rate[12][9] =  336.00;
    rate[12][10] =  527.00; rate[12][11] =  243.00; rate[13][0] =   18.00; rate[13][1] =   14.00;
    rate[13][2] =   14.00; rate[13][3] =    0.00; rate[13][4] =    0.00; rate[13][5] =    0.00;
    rate[13][6] =    0.00; rate[13][7] =   15.00; rate[13][8] =   48.00; rate[13][9] =  196.00;
    rate[13][10] =  157.00; rate[13][11] =    0.00; rate[13][12] =   92.00; rate[14][0] =  250.00;
    rate[14][1] =  103.00; rate[14][2] =   42.00; rate[14][3] =   13.00; rate[14][4] =   19.00;
    rate[14][5] =  153.00; rate[14][6] =   51.00; rate[14][7] =   34.00; rate[14][8] =   94.00;
    rate[14][9] =   12.00; rate[14][10] =   32.00; rate[14][11] =   33.00; rate[14][12] =   17.00;
    rate[14][13] =   11.00; rate[15][0] =  409.00; rate[15][1] =  154.00; rate[15][2] =  495.00;
    rate[15][3] =   95.00; rate[15][4] =  161.00; rate[15][5] =   56.00; rate[15][6] =   79.00;
    rate[15][7] =  234.00; rate[15][8] =   35.00; rate[15][9] =   24.00; rate[15][10] =   17.00;
    rate[15][11] =   96.00; rate[15][12] =   62.00; rate[15][13] =   46.00; rate[15][14] =  245.00;
    rate[16][0] =  371.00; rate[16][1] =   26.00; rate[16][2] =  229.00; rate[16][3] =   66.00;
    rate[16][4] =   16.00; rate[16][5] =   53.00; rate[16][6] =   34.00; rate[16][7] =   30.00;
    rate[16][8] =   22.00; rate[16][9] =  192.00; rate[16][10] =   33.00; rate[16][11] =  136.00;
    rate[16][12] =  104.00; rate[16][13] =   13.00; rate[16][14] =   78.00; rate[16][15] =  550.00;
    rate[17][0] =    0.00; rate[17][1] =  201.00; rate[17][2] =   23.00; rate[17][3] =    0.00;
    rate[17][4] =    0.00; rate[17][5] =    0.00; rate[17][6] =    0.00; rate[17][7] =    0.00;
    rate[17][8] =   27.00; rate[17][9] =    0.00; rate[17][10] =   46.00; rate[17][11] =    0.00;
    rate[17][12] =    0.00; rate[17][13] =   76.00; rate[17][14] =    0.00; rate[17][15] =   75.00;
    rate[17][16] =    0.00; rate[18][0] =   24.00; rate[18][1] =    8.00; rate[18][2] =   95.00;
    rate[18][3] =    0.00; rate[18][4] =   96.00; rate[18][5] =    0.00; rate[18][6] =   22.00;
    rate[18][7] =    0.00; rate[18][8] =  127.00; rate[18][9] =   37.00; rate[18][10] =   28.00;
    rate[18][11] =   13.00; rate[18][12] =    0.00; rate[18][13] =  698.00; rate[18][14] =    0.00;
    rate[18][15] =   34.00; rate[18][16] =   42.00; rate[18][17] =   61.00; rate[19][0] =  208.00;
    rate[19][1] =   24.00; rate[19][2] =   15.00; rate[19][3] =   18.00; rate[19][4] =   49.00;
    rate[19][5] =   35.00; rate[19][6] =   37.00; rate[19][7] =   54.00; rate[19][8] =   44.00;
    rate[19][9] =  889.00; rate[19][10] =  175.00; rate[19][11] =   10.00; rate[19][12] =  258.00;
    rate[19][13] =   12.00; rate[19][14] =   48.00; rate[19][15] =   30.00; rate[19][16] =  157.00;
    rate[19][17] =    0.00; rate[19][18] =   28.00;
    
       for (int i = 0; i < 20; i++) for (int j = i + 1; j < 20; j++) rate[i][j] = rate[j][i]; return rate;
    }

    @Override
    public double[] getEmpiricalFrequencies() {
        double[] f = new double[20];
    
    f[ 0] = 0.087127; f[ 1] = 0.040904; f[ 2] = 0.040432; f[ 3] = 0.046872;
    f[ 4] = 0.033474; f[ 5] = 0.038255; f[ 6] = 0.049530; f[ 7] = 0.088612;
    f[ 8] = 0.033618; f[ 9] = 0.036886; f[10] = 0.085357; f[11] = 0.080482;
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
