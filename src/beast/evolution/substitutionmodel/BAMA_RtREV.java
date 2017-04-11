package beast.evolution.substitutionmodel;

import beast.evolution.datatype.Aminoacid;
import beast.evolution.datatype.DataType;
/** model data from codonPHYML, which is based on PHYML **/

class BAMA_RtREV extends EmpiricalSubstitutionModel {
    @Override
    double[][] getEmpiricalRates() {
        double[][] rate = new double[20][20];

    /* 
     This model has been 'translated' from John Huelsenbeck and Fredrik Ronquist
     MrBayes program into PHYML format by Federico Abascal. Many thanks to them. 
     */
    
    /* 
     Dimmic M.W., J.S. Rest, D.P. Mindell, and D. Goldstein. 2002. RArtREV:
     An amino acid substitution matrix for inference of retrovirus and
     reverse transcriptase phyLOGeny. Journal of Molecular Evolution
     55: 65-73.
     */
    
    
    rate[1][0]= 34;         rate[2][0]= 51;         rate[2][1]= 35;         rate[3][0]= 10;         
    rate[3][1]= 30;         rate[3][2]= 384;        rate[4][0]= 439;        rate[4][1]= 92;         
    rate[4][2]= 128;        rate[4][3]= 1;          rate[5][0]= 32;         rate[5][1]= 221;        
    rate[5][2]= 236;        rate[5][3]= 78;         rate[5][4]= 70;         rate[6][0]= 81;         
    rate[6][1]= 10;         rate[6][2]= 79;         rate[6][3]= 542;        rate[6][4]= 1;          
    rate[6][5]= 372;        rate[7][0]= 135;        rate[7][1]= 41;         rate[7][2]= 94;         
    rate[7][3]= 61;         rate[7][4]= 48;         rate[7][5]= 18;         rate[7][6]= 70;         
    rate[8][0]= 30;         rate[8][1]= 90;         rate[8][2]= 320;        rate[8][3]= 91;         
    rate[8][4]= 124;        rate[8][5]= 387;        rate[8][6]= 34;         rate[8][7]= 68;         
    rate[9][0]= 1;          rate[9][1]= 24;         rate[9][2]= 35;         rate[9][3]= 1;          
    rate[9][4]= 104;        rate[9][5]= 33;         rate[9][6]= 1;          rate[9][7]= 1;          
    rate[9][8]= 34;         rate[10][0]= 45;        rate[10][1]= 18;        rate[10][2]= 15;        
    rate[10][3]= 5;         rate[10][4]= 110;       rate[10][5]= 54;        rate[10][6]= 21;        
    rate[10][7]= 3;         rate[10][8]= 51;        rate[10][9]= 385;       rate[11][0]= 38;        
    rate[11][1]= 593;       rate[11][2]= 123;       rate[11][3]= 20;        rate[11][4]= 16;        
    rate[11][5]= 309;       rate[11][6]= 141;       rate[11][7]= 30;        rate[11][8]= 76;        
    rate[11][9]= 34;        rate[11][10]= 23;       rate[12][0]= 235;       rate[12][1]= 57;        
    rate[12][2]= 1;         rate[12][3]= 1;         rate[12][4]= 156;       rate[12][5]= 158;       
    rate[12][6]= 1;         rate[12][7]= 37;        rate[12][8]= 116;       rate[12][9]= 375;       
    rate[12][10]= 581;      rate[12][11]= 134;      rate[13][0]= 1;         rate[13][1]= 7;         
    rate[13][2]= 49;        rate[13][3]= 1;         rate[13][4]= 70;        rate[13][5]= 1;         
    rate[13][6]= 1;         rate[13][7]= 7;         rate[13][8]= 141;       rate[13][9]= 64;        
    rate[13][10]= 179;      rate[13][11]= 14;       rate[13][12]= 247;      rate[14][0]= 97;        
    rate[14][1]= 24;        rate[14][2]= 33;        rate[14][3]= 55;        rate[14][4]= 1;         
    rate[14][5]= 68;        rate[14][6]= 52;        rate[14][7]= 17;        rate[14][8]= 44;        
    rate[14][9]= 10;        rate[14][10]= 22;       rate[14][11]= 43;       rate[14][12]= 1;        
    rate[14][13]= 11;       rate[15][0]= 460;       rate[15][1]= 102;       rate[15][2]= 294;       
    rate[15][3]= 136;       rate[15][4]= 75;        rate[15][5]= 225;       rate[15][6]= 95;        
    rate[15][7]= 152;       rate[15][8]= 183;       rate[15][9]= 4;         rate[15][10]= 24;       
    rate[15][11]= 77;       rate[15][12]= 1;        rate[15][13]= 20;       rate[15][14]= 134;      
    rate[16][0]= 258;       rate[16][1]= 64;        rate[16][2]= 148;       rate[16][3]= 55;        
    rate[16][4]= 117;       rate[16][5]= 146;       rate[16][6]= 82;        rate[16][7]= 7;         
    rate[16][8]= 49;        rate[16][9]= 72;        rate[16][10]= 25;       rate[16][11]= 110;      
    rate[16][12]= 131;      rate[16][13]= 69;       rate[16][14]= 62;       rate[16][15]= 671;      
    rate[17][0]= 5;         rate[17][1]= 13;        rate[17][2]= 16;        rate[17][3]= 1;         
    rate[17][4]= 55;        rate[17][5]= 10;        rate[17][6]= 17;        rate[17][7]= 23;        
    rate[17][8]= 48;        rate[17][9]= 39;        rate[17][10]= 47;       rate[17][11]= 6;        
    rate[17][12]= 111;      rate[17][13]= 182;      rate[17][14]= 9;        rate[17][15]= 14;       
    rate[17][16]= 1;        rate[18][0]= 55;        rate[18][1]= 47;        rate[18][2]= 28;        
    rate[18][3]= 1;         rate[18][4]= 131;       rate[18][5]= 45;        rate[18][6]= 1;         
    rate[18][7]= 21;        rate[18][8]= 307;       rate[18][9]= 26;        rate[18][10]= 64;       
    rate[18][11]= 1;        rate[18][12]= 74;       rate[18][13]= 1017;     rate[18][14]= 14;       
    rate[18][15]= 31;       rate[18][16]= 34;       rate[18][17]= 176;      rate[19][0]= 197;       
    rate[19][1]= 29;        rate[19][2]= 21;        rate[19][3]= 6;         rate[19][4]= 295;       
    rate[19][5]= 36;        rate[19][6]= 35;        rate[19][7]= 3;         rate[19][8]= 1;         
    rate[19][9]= 1048;      rate[19][10]= 112;      rate[19][11]= 19;       rate[19][12]= 236;      
    rate[19][13]= 92;       rate[19][14]= 25;       rate[19][15]= 39;       rate[19][16]= 196;      
    rate[19][17]= 26;       rate[19][18]= 59;       
    
       return rate;
    }

    @Override
    public double[] getEmpiricalFrequencies() {
        double[] f = new double[20];
    
    f[0]= 0.0646;           f[1]= 0.0453;           f[2]= 0.0376;           f[3]= 0.0422;           
    f[4]= 0.0114;           f[5]= 0.0606;           f[6]= 0.0607;           f[7]= 0.0639;           
    f[8]= 0.0273;           f[9]= 0.0679;           f[10]= 0.1018;          f[11]= 0.0751;          
    f[12]= 0.015;           f[13]= 0.0287;          f[14]= 0.0681;          f[15]= 0.0488;          
    f[16]= 0.0622;          f[17]= 0.0251;          f[18]= 0.0318;          f[19]= 0.0619;          
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
