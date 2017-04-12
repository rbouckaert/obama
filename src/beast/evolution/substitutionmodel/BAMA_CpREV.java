package beast.evolution.substitutionmodel;

import beast.evolution.datatype.Aminoacid;
import beast.evolution.datatype.DataType;
/** model data from codonPHYML, which is based on PHYML **/

public class BAMA_CpREV extends EmpiricalSubstitutionModel {
    @Override
    double[][] getEmpiricalRates() {
        double[][] rate = new double[20][20];

    /* 
     This model has been 'translated' from John Huelsenbeck and Fredrik Ronquist
     MrBayes program into PHYML format by Federico Abascal. Many thanks to them. 
     */
    
    /*
     Adachi, J., P. Waddell, W. Martin, and M. Hasegawa. 2000. Plastid          
     genome phyLOGeny and a model of amino acid substitution for proteins    
     encoded by chloroplast DNA. Journal of Molecular Evolution              
     50:348-358.
     */
    
    
    rate[1][0]= 105;        rate[2][0]= 227;        rate[2][1]= 357;        rate[3][0]= 175;        
    rate[3][1]= 43;         rate[3][2]= 4435;       rate[4][0]= 669;        rate[4][1]= 823;        
    rate[4][2]= 538;        rate[4][3]= 10;         rate[5][0]= 157;        rate[5][1]= 1745;       
    rate[5][2]= 768;        rate[5][3]= 400;        rate[5][4]= 10;         rate[6][0]= 499;        
    rate[6][1]= 152;        rate[6][2]= 1055;       rate[6][3]= 3691;       rate[6][4]= 10;         
    rate[6][5]= 3122;       rate[7][0]= 665;        rate[7][1]= 243;        rate[7][2]= 653;        
    rate[7][3]= 431;        rate[7][4]= 303;        rate[7][5]= 133;        rate[7][6]= 379;        
    rate[8][0]= 66;         rate[8][1]= 715;        rate[8][2]= 1405;       rate[8][3]= 331;        
    rate[8][4]= 441;        rate[8][5]= 1269;       rate[8][6]= 162;        rate[8][7]= 19;         
    rate[9][0]= 145;        rate[9][1]= 136;        rate[9][2]= 168;        rate[9][3]= 10;         
    rate[9][4]= 280;        rate[9][5]= 92;         rate[9][6]= 148;        rate[9][7]= 40;         
    rate[9][8]= 29;         rate[10][0]= 197;       rate[10][1]= 203;       rate[10][2]= 113;       
    rate[10][3]= 10;        rate[10][4]= 396;       rate[10][5]= 286;       rate[10][6]= 82;        
    rate[10][7]= 20;        rate[10][8]= 66;        rate[10][9]= 1745;      rate[11][0]= 236;       
    rate[11][1]= 4482;      rate[11][2]= 2430;      rate[11][3]= 412;       rate[11][4]= 48;        
    rate[11][5]= 3313;      rate[11][6]= 2629;      rate[11][7]= 263;       rate[11][8]= 305;       
    rate[11][9]= 345;       rate[11][10]= 218;      rate[12][0]= 185;       rate[12][1]= 125;       
    rate[12][2]= 61;        rate[12][3]= 47;        rate[12][4]= 159;       rate[12][5]= 202;       
    rate[12][6]= 113;       rate[12][7]= 21;        rate[12][8]= 10;        rate[12][9]= 1772;      
    rate[12][10]= 1351;     rate[12][11]= 193;      rate[13][0]= 68;        rate[13][1]= 53;        
    rate[13][2]= 97;        rate[13][3]= 22;        rate[13][4]= 726;       rate[13][5]= 10;        
    rate[13][6]= 145;       rate[13][7]= 25;        rate[13][8]= 127;       rate[13][9]= 454;       
    rate[13][10]= 1268;     rate[13][11]= 72;       rate[13][12]= 327;      rate[14][0]= 490;       
    rate[14][1]= 87;        rate[14][2]= 173;       rate[14][3]= 170;       rate[14][4]= 285;       
    rate[14][5]= 323;       rate[14][6]= 185;       rate[14][7]= 28;        rate[14][8]= 152;       
    rate[14][9]= 117;       rate[14][10]= 219;      rate[14][11]= 302;      rate[14][12]= 100;      
    rate[14][13]= 43;       rate[15][0]= 2440;      rate[15][1]= 385;       rate[15][2]= 2085;      
    rate[15][3]= 590;       rate[15][4]= 2331;      rate[15][5]= 396;       rate[15][6]= 568;       
    rate[15][7]= 691;       rate[15][8]= 303;       rate[15][9]= 216;       rate[15][10]= 516;      
    rate[15][11]= 868;      rate[15][12]= 93;       rate[15][13]= 487;      rate[15][14]= 1202;     
    rate[16][0]= 1340;      rate[16][1]= 314;       rate[16][2]= 1393;      rate[16][3]= 266;       
    rate[16][4]= 576;       rate[16][5]= 241;       rate[16][6]= 369;       rate[16][7]= 92;        
    rate[16][8]= 32;        rate[16][9]= 1040;      rate[16][10]= 156;      rate[16][11]= 918;      
    rate[16][12]= 645;      rate[16][13]= 148;      rate[16][14]= 260;      rate[16][15]= 2151;     
    rate[17][0]= 14;        rate[17][1]= 230;       rate[17][2]= 40;        rate[17][3]= 18;        
    rate[17][4]= 435;       rate[17][5]= 53;        rate[17][6]= 63;        rate[17][7]= 82;        
    rate[17][8]= 69;        rate[17][9]= 42;        rate[17][10]= 159;      rate[17][11]= 10;       
    rate[17][12]= 86;       rate[17][13]= 468;      rate[17][14]= 49;       rate[17][15]= 73;       
    rate[17][16]= 29;       rate[18][0]= 56;        rate[18][1]= 323;       rate[18][2]= 754;       
    rate[18][3]= 281;       rate[18][4]= 1466;      rate[18][5]= 391;       rate[18][6]= 142;       
    rate[18][7]= 10;        rate[18][8]= 1971;      rate[18][9]= 89;        rate[18][10]= 189;      
    rate[18][11]= 247;      rate[18][12]= 215;      rate[18][13]= 2370;     rate[18][14]= 97;       
    rate[18][15]= 522;      rate[18][16]= 71;       rate[18][17]= 346;      rate[19][0]= 968;       
    rate[19][1]= 92;        rate[19][2]= 83;        rate[19][3]= 75;        rate[19][4]= 592;       
    rate[19][5]= 54;        rate[19][6]= 200;       rate[19][7]= 91;        rate[19][8]= 25;        
    rate[19][9]= 4797;      rate[19][10]= 865;      rate[19][11]= 249;      rate[19][12]= 475;      
    rate[19][13]= 317;      rate[19][14]= 122;      rate[19][15]= 167;      rate[19][16]= 760;      
    rate[19][17]= 10;       rate[19][18]= 119;      
    
       for (int i = 0; i < 20; i++) for (int j = i + 1; j < 20; j++) rate[i][j] = rate[j][i]; return rate;
    }

    @Override
    public double[] getEmpiricalFrequencies() {
        double[] f = new double[20];
    
    f[0]= 0.076;            f[1]= 0.062;            f[2]= 0.041;            f[3]= 0.037;            
    f[4]= 0.009;            f[5]= 0.038;            f[6]= 0.049;            f[7]= 0.084;            
    f[8]= 0.025;            f[9]= 0.081;            f[10]= 0.101;           f[11]= 0.05;            
    f[12]= 0.022;           f[13]= 0.051;           f[14]= 0.043;           f[15]= 0.062;           
    f[16]= 0.054;           f[17]= 0.018;           f[18]= 0.031;           f[19]= 0.066;           
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
