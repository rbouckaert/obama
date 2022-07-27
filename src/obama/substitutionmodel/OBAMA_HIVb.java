package obama.substitutionmodel;

import beast.base.core.Description;
import beast.base.evolution.datatype.Aminoacid;
import beast.base.evolution.datatype.DataType;
import beast.base.evolution.substitutionmodel.EmpiricalSubstitutionModel;
/** model data from codonPHYML, which is based on PHYML **/

@Description("HIVb substitution model for amino acids")
public class OBAMA_HIVb extends EmpiricalSubstitutionModel {
    @Override
    public double[][] getEmpiricalRates() {
        double[][] rate = new double[20][20];

    /*
    //added by FEDE
     Nickle DC, Heath L, Jensen MA, Gilbert PB, Mullins JI, Kosakovsky Pond SL.
     HIV-Specific Probabilistic Models of Protein Evolution.
     PLoS ONE. 2007 Jun 6;2:e503.
     
     [thanks to Sergei L. Kosakovsky]
     
     Translated from HYPHY to Phyml format by Federico Abascal.
     */
    rate[1][0]= 0.307507;        rate[2][0]= 0.005;           rate[2][1]= 0.295543;        rate[3][0]= 1.45504;         
    rate[3][1]= 0.005;           rate[3][2]= 17.6612;         rate[4][0]= 0.123758;        rate[4][1]= 0.351721;        
    rate[4][2]= 0.0860642;       rate[4][3]= 0.005;           rate[5][0]= 0.0551128;       rate[5][1]= 3.4215;          
    rate[5][2]= 0.672052;        rate[5][3]= 0.005;           rate[5][4]= 0.005;           rate[6][0]= 1.48135;         
    rate[6][1]= 0.0749218;       rate[6][2]= 0.0792633;       rate[6][3]= 10.5872;         rate[6][4]= 0.005;           
    rate[6][5]= 2.5602;          rate[7][0]= 2.13536;         rate[7][1]= 3.65345;         rate[7][2]= 0.323401;        
    rate[7][3]= 2.83806;         rate[7][4]= 0.897871;        rate[7][5]= 0.0619137;       rate[7][6]= 3.92775;         
    rate[8][0]= 0.0847613;       rate[8][1]= 9.04044;         rate[8][2]= 7.64585;         rate[8][3]= 1.9169;          
    rate[8][4]= 0.240073;        rate[8][5]= 7.05545;         rate[8][6]= 0.11974;         rate[8][7]= 0.005;           
    rate[9][0]= 0.005;           rate[9][1]= 0.677289;        rate[9][2]= 0.680565;        rate[9][3]= 0.0176792;       
    rate[9][4]= 0.005;           rate[9][5]= 0.005;           rate[9][6]= 0.00609079;      rate[9][7]= 0.005;           
    rate[9][8]= 0.103111;        rate[10][0]= 0.215256;       rate[10][1]= 0.701427;       rate[10][2]= 0.005;          
    rate[10][3]= 0.00876048;     rate[10][4]= 0.129777;       rate[10][5]= 1.49456;        rate[10][6]= 0.005;          
    rate[10][7]= 0.005;          rate[10][8]= 1.74171;        rate[10][9]= 5.95879;        rate[11][0]= 0.005;          
    rate[11][1]= 20.45;          rate[11][2]= 7.90443;        rate[11][3]= 0.005;          rate[11][4]= 0.005;          
    rate[11][5]= 6.54737;        rate[11][6]= 4.61482;        rate[11][7]= 0.521705;       rate[11][8]= 0.005;          
    rate[11][9]= 0.322319;       rate[11][10]= 0.0814995;     rate[12][0]= 0.0186643;      rate[12][1]= 2.51394;        
    rate[12][2]= 0.005;          rate[12][3]= 0.005;          rate[12][4]= 0.005;          rate[12][5]= 0.303676;       
    rate[12][6]= 0.175789;       rate[12][7]= 0.005;          rate[12][8]= 0.005;          rate[12][9]= 11.2065;        
    rate[12][10]= 5.31961;       rate[12][11]= 1.28246;       rate[13][0]= 0.0141269;      rate[13][1]= 0.005;          
    rate[13][2]= 0.005;          rate[13][3]= 0.005;          rate[13][4]= 9.29815;        rate[13][5]= 0.005;          
    rate[13][6]= 0.005;          rate[13][7]= 0.291561;       rate[13][8]= 0.145558;       rate[13][9]= 3.39836;        
    rate[13][10]= 8.52484;       rate[13][11]= 0.0342658;     rate[13][12]= 0.188025;      rate[14][0]= 2.12217;        
    rate[14][1]= 1.28355;        rate[14][2]= 0.00739578;     rate[14][3]= 0.0342658;      rate[14][4]= 0.005;          
    rate[14][5]= 4.47211;        rate[14][6]= 0.0120226;      rate[14][7]= 0.005;          rate[14][8]= 2.45318;        
    rate[14][9]= 0.0410593;      rate[14][10]= 2.07757;       rate[14][11]= 0.0313862;     rate[14][12]= 0.005;         
    rate[14][13]= 0.005;         rate[15][0]= 2.46633;        rate[15][1]= 3.4791;         rate[15][2]= 13.1447;        
    rate[15][3]= 0.52823;        rate[15][4]= 4.69314;        rate[15][5]= 0.116311;       rate[15][6]= 0.005;          
    rate[15][7]= 4.38041;        rate[15][8]= 0.382747;       rate[15][9]= 1.21803;        rate[15][10]= 0.927656;      
    rate[15][11]= 0.504111;      rate[15][12]= 0.005;         rate[15][13]= 0.956472;      rate[15][14]= 5.37762;       
    rate[16][0]= 15.9183;        rate[16][1]= 2.86868;        rate[16][2]= 6.88667;        rate[16][3]= 0.274724;       
    rate[16][4]= 0.739969;       rate[16][5]= 0.243589;       rate[16][6]= 0.289774;       rate[16][7]= 0.369615;       
    rate[16][8]= 0.711594;       rate[16][9]= 8.61217;        rate[16][10]= 0.0437673;     rate[16][11]= 4.67142;       
    rate[16][12]= 4.94026;       rate[16][13]= 0.0141269;     rate[16][14]= 2.01417;       rate[16][15]= 8.93107;       
    rate[17][0]= 0.005;          rate[17][1]= 0.991338;       rate[17][2]= 0.005;          rate[17][3]= 0.005;          
    rate[17][4]= 2.63277;        rate[17][5]= 0.026656;       rate[17][6]= 0.005;          rate[17][7]= 1.21674;        
    rate[17][8]= 0.0695179;      rate[17][9]= 0.005;          rate[17][10]= 0.748843;      rate[17][11]= 0.005;         
    rate[17][12]= 0.089078;      rate[17][13]= 0.829343;      rate[17][14]= 0.0444506;     rate[17][15]= 0.0248728;     
    rate[17][16]= 0.005;         rate[18][0]= 0.005;          rate[18][1]= 0.00991826;     rate[18][2]= 1.76417;        
    rate[18][3]= 0.674653;       rate[18][4]= 7.57932;        rate[18][5]= 0.113033;       rate[18][6]= 0.0792633;      
    rate[18][7]= 0.005;          rate[18][8]= 18.6943;        rate[18][9]= 0.148168;       rate[18][10]= 0.111986;      
    rate[18][11]= 0.005;         rate[18][12]= 0.005;         rate[18][13]= 15.34;         rate[18][14]= 0.0304381;     
    rate[18][15]= 0.648024;      rate[18][16]= 0.105652;      rate[18][17]= 1.28022;       rate[19][0]= 7.61428;        
    rate[19][1]= 0.0812454;      rate[19][2]= 0.026656;       rate[19][3]= 1.04793;        rate[19][4]= 0.420027;       
    rate[19][5]= 0.0209153;      rate[19][6]= 1.02847;        rate[19][7]= 0.953155;       rate[19][8]= 0.005;          
    rate[19][9]= 17.7389;        rate[19][10]= 1.41036;       rate[19][11]= 0.265829;      rate[19][12]= 6.8532;        
    rate[19][13]= 0.723274;      rate[19][14]= 0.005;         rate[19][15]= 0.0749218;     rate[19][16]= 0.709226;      
    rate[19][17]= 0.005;         rate[19][18]= 0.0410593;     
       for (int i = 0; i < 20; i++) for (int j = i + 1; j < 20; j++) rate[i][j] = rate[j][i]; return rate;
    }

    @Override
    public double[] getEmpiricalFrequencies() {
        double[] f = new double[20];
    
    f[0]= 0.060490222;           f[1]= 0.066039665;           f[2]= 0.044127815;           f[3]= 0.042109048;           
    f[4]= 0.020075899;           f[5]= 0.053606488;           f[6]= 0.071567447;           f[7]= 0.072308239;           
    f[8]= 0.022293943;           f[9]= 0.069730629;           f[10]= 0.098851122;          f[11]= 0.056968211;          
    f[12]= 0.019768318;          f[13]= 0.028809447;          f[14]= 0.046025282;          f[15]= 0.05060433;           
    f[16]= 0.053636813;          f[17]= 0.033011601;          f[18]= 0.028350243;          f[19]= 0.061625237;          
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
