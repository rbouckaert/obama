package beast.evolution.substitutionmodel;

import beast.evolution.datatype.Aminoacid;
import beast.evolution.datatype.DataType;
/** model data from codonPHYML, which is based on PHYML **/

public class OBAMA_HIVw extends EmpiricalSubstitutionModel {
    @Override
    double[][] getEmpiricalRates() {
        double[][] rate = new double[20][20];

    /*
     Nickle DC, Heath L, Jensen MA, Gilbert PB, Mullins JI, Kosakovsky Pond SL.
     HIV-Specific Probabilistic Models of Protein Evolution.
     PLoS ONE. 2007 Jun 6;2:e503.
     
     [thanks to Sergei L. Kosakovsky]
     
     Translated from HYPHY to Phyml format by Federico Abascal.
     */
    
    
    rate[1][0]= 0.0744808;       rate[2][0]= 0.617509;        rate[2][1]= 0.16024;         rate[3][0]= 4.43521;         
    rate[3][1]= 0.0674539;       rate[3][2]= 29.4087;         rate[4][0]= 0.167653;        rate[4][1]= 2.86364;         
    rate[4][2]= 0.0604932;       rate[4][3]= 0.005;           rate[5][0]= 0.005;           rate[5][1]= 10.6746;         
    rate[5][2]= 0.342068;        rate[5][3]= 0.005;           rate[5][4]= 0.005;           rate[6][0]= 5.56325;         
    rate[6][1]= 0.0251632;       rate[6][2]= 0.201526;        rate[6][3]= 12.1233;         rate[6][4]= 0.005;           
    rate[6][5]= 3.20656;         rate[7][0]= 1.8685;          rate[7][1]= 13.4379;         rate[7][2]= 0.0604932;       
    rate[7][3]= 10.3969;         rate[7][4]= 0.0489798;       rate[7][5]= 0.0604932;       rate[7][6]= 14.7801;         
    rate[8][0]= 0.005;           rate[8][1]= 6.84405;         rate[8][2]= 8.59876;         rate[8][3]= 2.31779;         
    rate[8][4]= 0.005;           rate[8][5]= 18.5465;         rate[8][6]= 0.005;           rate[8][7]= 0.005;           
    rate[9][0]= 0.005;           rate[9][1]= 1.34069;         rate[9][2]= 0.987028;        rate[9][3]= 0.145124;        
    rate[9][4]= 0.005;           rate[9][5]= 0.0342252;       rate[9][6]= 0.0390512;       rate[9][7]= 0.005;           
    rate[9][8]= 0.005;           rate[10][0]= 0.16024;        rate[10][1]= 0.586757;       rate[10][2]= 0.005;          
    rate[10][3]= 0.005;          rate[10][4]= 0.005;          rate[10][5]= 2.89048;        rate[10][6]= 0.129839;       
    rate[10][7]= 0.0489798;      rate[10][8]= 1.76382;        rate[10][9]= 9.10246;        rate[11][0]= 0.592784;       
    rate[11][1]= 39.8897;        rate[11][2]= 10.6655;        rate[11][3]= 0.894313;       rate[11][4]= 0.005;          
    rate[11][5]= 13.0705;        rate[11][6]= 23.9626;        rate[11][7]= 0.279425;       rate[11][8]= 0.22406;        
    rate[11][9]= 0.817481;       rate[11][10]= 0.005;         rate[12][0]= 0.005;          rate[12][1]= 3.28652;        
    rate[12][2]= 0.201526;       rate[12][3]= 0.005;          rate[12][4]= 0.005;          rate[12][5]= 0.005;          
    rate[12][6]= 0.005;          rate[12][7]= 0.0489798;      rate[12][8]= 0.005;          rate[12][9]= 17.3064;        
    rate[12][10]= 11.3839;       rate[12][11]= 4.09564;       rate[13][0]= 0.597923;       rate[13][1]= 0.005;          
    rate[13][2]= 0.005;          rate[13][3]= 0.005;          rate[13][4]= 0.362959;       rate[13][5]= 0.005;          
    rate[13][6]= 0.005;          rate[13][7]= 0.005;          rate[13][8]= 0.005;          rate[13][9]= 1.48288;        
    rate[13][10]= 7.48781;       rate[13][11]= 0.005;         rate[13][12]= 0.005;         rate[14][0]= 1.00981;        
    rate[14][1]= 0.404723;       rate[14][2]= 0.344848;       rate[14][3]= 0.005;          rate[14][4]= 0.005;          
    rate[14][5]= 3.04502;        rate[14][6]= 0.005;          rate[14][7]= 0.005;          rate[14][8]= 13.9444;        
    rate[14][9]= 0.005;          rate[14][10]= 9.83095;       rate[14][11]= 0.111928;      rate[14][12]= 0.005;         
    rate[14][13]= 0.0342252;     rate[15][0]= 8.5942;         rate[15][1]= 8.35024;        rate[15][2]= 14.5699;        
    rate[15][3]= 0.427881;       rate[15][4]= 1.12195;        rate[15][5]= 0.16024;        rate[15][6]= 0.005;          
    rate[15][7]= 6.27966;        rate[15][8]= 0.725157;       rate[15][9]= 0.740091;       rate[15][10]= 6.14396;       
    rate[15][11]= 0.005;         rate[15][12]= 0.392575;      rate[15][13]= 4.27939;       rate[15][14]= 14.249;        
    rate[16][0]= 24.1422;        rate[16][1]= 0.928203;       rate[16][2]= 4.54206;        rate[16][3]= 0.630395;       
    rate[16][4]= 0.005;          rate[16][5]= 0.203091;       rate[16][6]= 0.458743;       rate[16][7]= 0.0489798;      
    rate[16][8]= 0.95956;        rate[16][9]= 9.36345;        rate[16][10]= 0.005;         rate[16][11]= 4.04802;       
    rate[16][12]= 7.41313;       rate[16][13]= 0.114512;      rate[16][14]= 4.33701;       rate[16][15]= 6.34079;       
    rate[17][0]= 0.005;          rate[17][1]= 5.96564;        rate[17][2]= 0.005;          rate[17][3]= 0.005;          
    rate[17][4]= 5.49894;        rate[17][5]= 0.0443298;      rate[17][6]= 0.005;          rate[17][7]= 2.8258;         
    rate[17][8]= 0.005;          rate[17][9]= 0.005;          rate[17][10]= 1.37031;       rate[17][11]= 0.005;         
    rate[17][12]= 0.005;         rate[17][13]= 0.005;         rate[17][14]= 0.005;         rate[17][15]= 1.10156;       
    rate[17][16]= 0.005;         rate[18][0]= 0.005;          rate[18][1]= 0.005;          rate[18][2]= 5.06475;        
    rate[18][3]= 2.28154;        rate[18][4]= 8.34835;        rate[18][5]= 0.005;          rate[18][6]= 0.005;          
    rate[18][7]= 0.005;          rate[18][8]= 47.4889;        rate[18][9]= 0.114512;       rate[18][10]= 0.005;         
    rate[18][11]= 0.005;         rate[18][12]= 0.579198;      rate[18][13]= 4.12728;       rate[18][14]= 0.005;         
    rate[18][15]= 0.933142;      rate[18][16]= 0.490608;      rate[18][17]= 0.005;         rate[19][0]= 24.8094;        
    rate[19][1]= 0.279425;       rate[19][2]= 0.0744808;      rate[19][3]= 2.91786;        rate[19][4]= 0.005;          
    rate[19][5]= 0.005;          rate[19][6]= 2.19952;        rate[19][7]= 2.79622;        rate[19][8]= 0.827479;       
    rate[19][9]= 24.8231;        rate[19][10]= 2.95344;       rate[19][11]= 0.128065;      rate[19][12]= 14.7683;       
    rate[19][13]= 2.28;          rate[19][14]= 0.005;         rate[19][15]= 0.862637;      rate[19][16]= 0.005;         
    rate[19][17]= 0.005;         rate[19][18]= 1.35482;       
       for (int i = 0; i < 20; i++) for (int j = i + 1; j < 20; j++) rate[i][j] = rate[j][i]; return rate;
    }

    @Override
    public double[] getEmpiricalFrequencies() {
        double[] f = new double[20];
    
    f[0]= 0.0377494;             f[1]= 0.057321;              f[2]= 0.0891129;             f[3]= 0.0342034;             
    f[4]= 0.0240105;             f[5]= 0.0437824;             f[6]= 0.0618606;             f[7]= 0.0838496;             
    f[8]= 0.0156076;             f[9]= 0.0983641;             f[10]= 0.0577867;            f[11]= 0.0641682;            
    f[12]= 0.0158419;            f[13]= 0.0422741;            f[14]= 0.0458601;            f[15]= 0.0550846;            
    f[16]= 0.0813774;            f[17]= 0.019597;             f[18]= 0.0205847;            f[19]= 0.0515639;            
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
