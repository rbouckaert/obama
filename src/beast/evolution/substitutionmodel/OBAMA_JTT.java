package beast.evolution.substitutionmodel;

import beast.base.core.Description;
import beast.base.evolution.datatype.Aminoacid;
import beast.base.evolution.datatype.DataType;
import beast.base.evolution.substitutionmodel.EmpiricalSubstitutionModel;
/** model data from codonPHYML, which is based on PHYML **/

@Description("JTT substitution model for amino acids")
public class OBAMA_JTT extends EmpiricalSubstitutionModel {
    @Override
    public double[][] getEmpiricalRates() {
        double[][] rate = new double[20][20];

    
    /* JTT's model data
     * D.T.Jones, W.R.Taylor and J.M.Thornton
     * "The rapid generation of mutation data matrices from protein sequences"
     * CABIOS  vol.8 no.3 1992 pp275-282
     */
    

    
    rate[1][0] =   58.00; rate[2][0] =   54.00; rate[2][1] =   45.00; rate[3][0] =   81.00;
    rate[3][1] =   16.00; rate[3][2] =  528.00; rate[4][0] =   56.00; rate[4][1] =  113.00;
    rate[4][2] =   34.00; rate[4][3] =   10.00; rate[5][0] =   57.00; rate[5][1] =  310.00;
    rate[5][2] =   86.00; rate[5][3] =   49.00; rate[5][4] =    9.00; rate[6][0] =  105.00;
    rate[6][1] =   29.00; rate[6][2] =   58.00; rate[6][3] =  767.00; rate[6][4] =    5.00;
    rate[6][5] =  323.00; rate[7][0] =  179.00; rate[7][1] =  137.00; rate[7][2] =   81.00;
    rate[7][3] =  130.00; rate[7][4] =   59.00; rate[7][5] =   26.00; rate[7][6] =  119.00;
    rate[8][0] =   27.00; rate[8][1] =  328.00; rate[8][2] =  391.00; rate[8][3] =  112.00;
    rate[8][4] =   69.00; rate[8][5] =  597.00; rate[8][6] =   26.00; rate[8][7] =   23.00;
    rate[9][0] =   36.00; rate[9][1] =   22.00; rate[9][2] =   47.00; rate[9][3] =   11.00;
    rate[9][4] =   17.00; rate[9][5] =    9.00; rate[9][6] =   12.00; rate[9][7] =    6.00;
    rate[9][8] =   16.00; rate[10][0] =   30.00; rate[10][1] =   38.00; rate[10][2] =   12.00;
    rate[10][3] =    7.00; rate[10][4] =   23.00; rate[10][5] =   72.00; rate[10][6] =    9.00;
    rate[10][7] =    6.00; rate[10][8] =   56.00; rate[10][9] =  229.00; rate[11][0] =   35.00;
    rate[11][1] =  646.00; rate[11][2] =  263.00; rate[11][3] =   26.00; rate[11][4] =    7.00;
    rate[11][5] =  292.00; rate[11][6] =  181.00; rate[11][7] =   27.00; rate[11][8] =   45.00;
    rate[11][9] =   21.00; rate[11][10] =   14.00; rate[12][0] =   54.00; rate[12][1] =   44.00;
    rate[12][2] =   30.00; rate[12][3] =   15.00; rate[12][4] =   31.00; rate[12][5] =   43.00;
    rate[12][6] =   18.00; rate[12][7] =   14.00; rate[12][8] =   33.00; rate[12][9] =  479.00;
    rate[12][10] =  388.00; rate[12][11] =   65.00; rate[13][0] =   15.00; rate[13][1] =    5.00;
    rate[13][2] =   10.00; rate[13][3] =    4.00; rate[13][4] =   78.00; rate[13][5] =    4.00;
    rate[13][6] =    5.00; rate[13][7] =    5.00; rate[13][8] =   40.00; rate[13][9] =   89.00;
    rate[13][10] =  248.00; rate[13][11] =    4.00; rate[13][12] =   43.00; rate[14][0] =  194.00;
    rate[14][1] =   74.00; rate[14][2] =   15.00; rate[14][3] =   15.00; rate[14][4] =   14.00;
    rate[14][5] =  164.00; rate[14][6] =   18.00; rate[14][7] =   24.00; rate[14][8] =  115.00;
    rate[14][9] =   10.00; rate[14][10] =  102.00; rate[14][11] =   21.00; rate[14][12] =   16.00;
    rate[14][13] =   17.00; rate[15][0] =  378.00; rate[15][1] =  101.00; rate[15][2] =  503.00;
    rate[15][3] =   59.00; rate[15][4] =  223.00; rate[15][5] =   53.00; rate[15][6] =   30.00;
    rate[15][7] =  201.00; rate[15][8] =   73.00; rate[15][9] =   40.00; rate[15][10] =   59.00;
    rate[15][11] =   47.00; rate[15][12] =   29.00; rate[15][13] =   92.00; rate[15][14] =  285.00;
    rate[16][0] =  475.00; rate[16][1] =   64.00; rate[16][2] =  232.00; rate[16][3] =   38.00;
    rate[16][4] =   42.00; rate[16][5] =   51.00; rate[16][6] =   32.00; rate[16][7] =   33.00;
    rate[16][8] =   46.00; rate[16][9] =  245.00; rate[16][10] =   25.00; rate[16][11] =  103.00;
    rate[16][12] =  226.00; rate[16][13] =   12.00; rate[16][14] =  118.00; rate[16][15] =  477.00;
    rate[17][0] =    9.00; rate[17][1] =  126.00; rate[17][2] =    8.00; rate[17][3] =    4.00;
    rate[17][4] =  115.00; rate[17][5] =   18.00; rate[17][6] =   10.00; rate[17][7] =   55.00;
    rate[17][8] =    8.00; rate[17][9] =    9.00; rate[17][10] =   52.00; rate[17][11] =   10.00;
    rate[17][12] =   24.00; rate[17][13] =   53.00; rate[17][14] =    6.00; rate[17][15] =   35.00;
    rate[17][16] =   12.00; rate[18][0] =   11.00; rate[18][1] =   20.00; rate[18][2] =   70.00;
    rate[18][3] =   46.00; rate[18][4] =  209.00; rate[18][5] =   24.00; rate[18][6] =    7.00;
    rate[18][7] =    8.00; rate[18][8] =  573.00; rate[18][9] =   32.00; rate[18][10] =   24.00;
    rate[18][11] =    8.00; rate[18][12] =   18.00; rate[18][13] =  536.00; rate[18][14] =   10.00;
    rate[18][15] =   63.00; rate[18][16] =   21.00; rate[18][17] =   71.00; rate[19][0] =  298.00;
    rate[19][1] =   17.00; rate[19][2] =   16.00; rate[19][3] =   31.00; rate[19][4] =   62.00;
    rate[19][5] =   20.00; rate[19][6] =   45.00; rate[19][7] =   47.00; rate[19][8] =   11.00;
    rate[19][9] =  961.00; rate[19][10] =  180.00; rate[19][11] =   14.00; rate[19][12] =  323.00;
    rate[19][13] =   62.00; rate[19][14] =   23.00; rate[19][15] =   38.00; rate[19][16] =  112.00;
    rate[19][17] =   25.00; rate[19][18] =   16.00;
    
       for (int i = 0; i < 20; i++) for (int j = i + 1; j < 20; j++) rate[i][j] = rate[j][i]; return rate;
    }

    @Override
    public double[] getEmpiricalFrequencies() {
        double[] f = new double[20];
    
    f[ 0] = 0.076748; f[ 1] = 0.051691; f[ 2] = 0.042645; f[ 3] = 0.051544;
    f[ 4] = 0.019803; f[ 5] = 0.040752; f[ 6] = 0.061830; f[ 7] = 0.073152;
    f[ 8] = 0.022944; f[ 9] = 0.053761; f[10] = 0.091904; f[11] = 0.058676;
    f[12] = 0.023826; f[13] = 0.040126; f[14] = 0.050901; f[15] = 0.068765;
    f[16] = 0.058565; f[17] = 0.014261; f[18] = 0.032102; f[19] = 0.066004; // was 0.066005, but does not add to 1
    
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
