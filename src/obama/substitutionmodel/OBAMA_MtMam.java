package obama.substitutionmodel;

import beast.base.core.Description;
import beast.base.evolution.datatype.Aminoacid;
import beast.base.evolution.datatype.DataType;
import beast.base.evolution.substitutionmodel.EmpiricalSubstitutionModel;
/** model data from codonPHYML, which is based on PHYML **/

@Description("MtMam substitution model for amino acids")
public class OBAMA_MtMam extends EmpiricalSubstitutionModel {
    @Override
    public double[][] getEmpiricalRates() {
        double[][] rate = new double[20][20];

    /* 
     This model has been 'translated' from Ziheng Yang's PAML program 
     into PHYML format by Federico Abascal. Many thanks to them. 
     */
    
    /*
     Cao, Y. et al. 1998 Conflict amongst individual mitochondrial 
     proteins in resolving the phyLOGeny of eutherian orders. Journal 
     of Molecular Evolution 15:1600-1611.
     */
    
    
    rate[1][0]= 32;              rate[2][0]= 2;    rate[2][1]= 4;               rate[3][0]= 11;
    rate[3][1]= 0;               rate[3][2]= 864;  rate[4][0]= 0;               rate[4][1]= 186;
    rate[4][2]= 0;               rate[4][3]= 0;    rate[5][0]= 0;               rate[5][1]= 246;
    rate[5][2]= 8;               rate[5][3]= 49;   rate[5][4]= 0;               rate[6][0]= 0;
    rate[6][1]= 0;               rate[6][2]= 0;    rate[6][3]= 569;             rate[6][4]= 0;
    rate[6][5]= 274;             rate[7][0]= 78;   rate[7][1]= 18;              rate[7][2]= 47;
    rate[7][3]= 79;              rate[7][4]= 0;    rate[7][5]= 0;               rate[7][6]= 22;
    rate[8][0]= 8;               rate[8][1]= 232;  rate[8][2]= 458;             rate[8][3]= 11;
    rate[8][4]= 305;             rate[8][5]= 550;  rate[8][6]= 22;              rate[8][7]= 0;
    rate[9][0]= 75;              rate[9][1]= 0;    rate[9][2]= 19;              rate[9][3]= 0;
    rate[9][4]= 41;              rate[9][5]= 0;    rate[9][6]= 0;               rate[9][7]= 0;
    rate[9][8]= 0;               rate[10][0]= 21;  rate[10][1]= 6;              rate[10][2]= 0;
    rate[10][3]= 0;              rate[10][4]= 27;  rate[10][5]= 20;             rate[10][6]= 0;
    rate[10][7]= 0;              rate[10][8]= 26;  rate[10][9]= 232;            rate[11][0]= 0;
    rate[11][1]= 50;             rate[11][2]= 408; rate[11][3]= 0;              rate[11][4]= 0;
    rate[11][5]= 242;            rate[11][6]= 215; rate[11][7]= 0;              rate[11][8]= 0;
    rate[11][9]= 6;              rate[11][10]= 4;  rate[12][0]= 76;             rate[12][1]= 0;
    rate[12][2]= 21;             rate[12][3]= 0;   rate[12][4]= 0;              rate[12][5]= 22;
    rate[12][6]= 0;              rate[12][7]= 0;   rate[12][8]= 0;              rate[12][9]= 378;
    rate[12][10]= 609;           rate[12][11]= 59; rate[13][0]= 0;              rate[13][1]= 0;
    rate[13][2]= 6;              rate[13][3]= 5;   rate[13][4]= 7;              rate[13][5]= 0;
    rate[13][6]= 0;              rate[13][7]= 0;   rate[13][8]= 0;              rate[13][9]= 57;
    rate[13][10]= 246;           rate[13][11]= 0;  rate[13][12]= 11;            rate[14][0]= 53;
    rate[14][1]= 9;              rate[14][2]= 33;  rate[14][3]= 2;              rate[14][4]= 0;
    rate[14][5]= 51;             rate[14][6]= 0;   rate[14][7]= 0;              rate[14][8]= 53;
    rate[14][9]= 5;              rate[14][10]= 43; rate[14][11]= 18;            rate[14][12]= 0;
    rate[14][13]= 17;            rate[15][0]= 342; rate[15][1]= 3;              rate[15][2]= 446;
    rate[15][3]= 16;             rate[15][4]= 347; rate[15][5]= 30;             rate[15][6]= 21;
    rate[15][7]= 112;            rate[15][8]= 20;  rate[15][9]= 0;              rate[15][10]= 74;
    rate[15][11]= 65;            rate[15][12]= 47; rate[15][13]= 90;            rate[15][14]= 202;
    rate[16][0]= 681;            rate[16][1]= 0;   rate[16][2]= 110;            rate[16][3]= 0;
    rate[16][4]= 114;            rate[16][5]= 0;   rate[16][6]= 4;              rate[16][7]= 0;
    rate[16][8]= 1;              rate[16][9]= 360; rate[16][10]= 34;            rate[16][11]= 50;
    rate[16][12]= 691;           rate[16][13]= 8;  rate[16][14]= 78;            rate[16][15]= 614;
    rate[17][0]= 5;              rate[17][1]= 16;  rate[17][2]= 6;              rate[17][3]= 0;
    rate[17][4]= 65;             rate[17][5]= 0;   rate[17][6]= 0;              rate[17][7]= 0;
    rate[17][8]= 0;              rate[17][9]= 0;   rate[17][10]= 12;            rate[17][11]= 0;
    rate[17][12]= 13;            rate[17][13]= 0;  rate[17][14]= 7;             rate[17][15]= 17;
    rate[17][16]= 0;             rate[18][0]= 0;   rate[18][1]= 0;              rate[18][2]= 156;
    rate[18][3]= 0;              rate[18][4]= 530; rate[18][5]= 54;             rate[18][6]= 0;
    rate[18][7]= 1;              rate[18][8]= 1525;rate[18][9]= 16;             rate[18][10]= 25;
    rate[18][11]= 67;            rate[18][12]= 0;  rate[18][13]= 682;           rate[18][14]= 8;
    rate[18][15]= 107;           rate[18][16]= 0;  rate[18][17]= 14;            rate[19][0]= 398;
    rate[19][1]= 0;              rate[19][2]= 0;   rate[19][3]= 10;             rate[19][4]= 0;
    rate[19][5]= 33;             rate[19][6]= 20;  rate[19][7]= 5;              rate[19][8]= 0;
    rate[19][9]= 2220;           rate[19][10]= 100;rate[19][11]= 0;             rate[19][12]= 832;
    rate[19][13]= 6;             rate[19][14]= 0;  rate[19][15]= 0;             rate[19][16]= 237;
    rate[19][17]= 0;             rate[19][18]= 0;
    
       for (int i = 0; i < 20; i++) for (int j = i + 1; j < 20; j++) rate[i][j] = rate[j][i]; return rate;
    }

    @Override
    public double[] getEmpiricalFrequencies() {
        double[] f = new double[20];
    
    f[0]= 0.0692;  f[1]=  0.0184;  f[2]= 0.04;    f[3]= 0.0186;
    f[4]= 0.0065;  f[5]=  0.0238;  f[6]= 0.0236;  f[7]= 0.0557;
    f[8]= 0.0277;  f[9]=  0.0905;  f[10]=0.1675;  f[11]= 0.0221;
    f[12]=0.0561;  f[13]= 0.0611;  f[14]=0.0536;  f[15]= 0.0725;
    f[16]=0.087;   f[17]= 0.0293;  f[18]=0.034;   f[19]= 0.0428;
    
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
