package beast.evolution.substitutionmodel;

import beast.evolution.datatype.Aminoacid;
import beast.evolution.datatype.DataType;
/** model data from codonPHYML, which is based on PHYML **/

class BAMA_VT extends EmpiricalSubstitutionModel {
    @Override
    double[][] getEmpiricalRates() {
        double[][] rate = new double[20][20];

    /* 
     This model has been 'translated' from John Huelsenbeck and Fredrik Ronquist
     MrBayes program into PHYML format by Federico Abascal. Many thanks to them. 
     */
    
    /*
     Muller, T., and M. Vingron. 2000. Modeling amino acid replacement.         
     Journal of Computational BioLOGy 7:761-776.                             
     */
    
    
    rate[1][0]= 0.233108;   rate[2][0]= 0.199097;   rate[2][1]= 0.210797;   rate[3][0]= 0.265145;   
    rate[3][1]= 0.105191;   rate[3][2]= 0.883422;   rate[4][0]= 0.227333;   rate[4][1]= 0.031726;   
    rate[4][2]= 0.027495;   rate[4][3]= 0.010313;   rate[5][0]= 0.310084;   rate[5][1]= 0.493763;   
    rate[5][2]= 0.2757;     rate[5][3]= 0.205842;   rate[5][4]= 0.004315;   rate[6][0]= 0.567957;   
    rate[6][1]= 0.25524;    rate[6][2]= 0.270417;   rate[6][3]= 1.599461;   rate[6][4]= 0.005321;   
    rate[6][5]= 0.960976;   rate[7][0]= 0.876213;   rate[7][1]= 0.156945;   rate[7][2]= 0.362028;   
    rate[7][3]= 0.311718;   rate[7][4]= 0.050876;   rate[7][5]= 0.12866;    rate[7][6]= 0.250447;   
    rate[8][0]= 0.078692;   rate[8][1]= 0.213164;   rate[8][2]= 0.290006;   rate[8][3]= 0.134252;   
    rate[8][4]= 0.016695;   rate[8][5]= 0.315521;   rate[8][6]= 0.104458;   rate[8][7]= 0.058131;   
    rate[9][0]= 0.222972;   rate[9][1]= 0.08151;    rate[9][2]= 0.087225;   rate[9][3]= 0.01172;    
    rate[9][4]= 0.046398;   rate[9][5]= 0.054602;   rate[9][6]= 0.046589;   rate[9][7]= 0.051089;   
    rate[9][8]= 0.020039;   rate[10][0]= 0.42463;   rate[10][1]= 0.192364;  rate[10][2]= 0.069245;  
    rate[10][3]= 0.060863;  rate[10][4]= 0.091709;  rate[10][5]= 0.24353;   rate[10][6]= 0.151924;  
    rate[10][7]= 0.087056;  rate[10][8]= 0.103552;  rate[10][9]= 2.08989;   rate[11][0]= 0.393245;  
    rate[11][1]= 1.755838;  rate[11][2]= 0.50306;   rate[11][3]= 0.261101;  rate[11][4]= 0.004067;  
    rate[11][5]= 0.738208;  rate[11][6]= 0.88863;   rate[11][7]= 0.193243;  rate[11][8]= 0.153323;  
    rate[11][9]= 0.093181;  rate[11][10]= 0.201204; rate[12][0]= 0.21155;   rate[12][1]= 0.08793;   
    rate[12][2]= 0.05742;   rate[12][3]= 0.012182;  rate[12][4]= 0.02369;   rate[12][5]= 0.120801;  
    rate[12][6]= 0.058643;  rate[12][7]= 0.04656;   rate[12][8]= 0.021157;  rate[12][9]= 0.493845;  
    rate[12][10]= 1.105667; rate[12][11]= 0.096474; rate[13][0]= 0.116646;  rate[13][1]= 0.042569;  
    rate[13][2]= 0.039769;  rate[13][3]= 0.016577;  rate[13][4]= 0.051127;  rate[13][5]= 0.026235;  
    rate[13][6]= 0.028168;  rate[13][7]= 0.050143;  rate[13][8]= 0.079807;  rate[13][9]= 0.32102;   
    rate[13][10]= 0.946499; rate[13][11]= 0.038261; rate[13][12]= 0.173052; rate[14][0]= 0.399143;  
    rate[14][1]= 0.12848;   rate[14][2]= 0.083956;  rate[14][3]= 0.160063;  rate[14][4]= 0.011137;  
    rate[14][5]= 0.15657;   rate[14][6]= 0.205134;  rate[14][7]= 0.124492;  rate[14][8]= 0.078892;  
    rate[14][9]= 0.054797;  rate[14][10]= 0.169784; rate[14][11]= 0.212302; rate[14][12]= 0.010363; 
    rate[14][13]= 0.042564; rate[15][0]= 1.817198;  rate[15][1]= 0.292327;  rate[15][2]= 0.847049;  
    rate[15][3]= 0.461519;  rate[15][4]= 0.17527;   rate[15][5]= 0.358017;  rate[15][6]= 0.406035;  
    rate[15][7]= 0.612843;  rate[15][8]= 0.167406;  rate[15][9]= 0.081567;  rate[15][10]= 0.214977; 
    rate[15][11]= 0.400072; rate[15][12]= 0.090515; rate[15][13]= 0.138119; rate[15][14]= 0.430431; 
    rate[16][0]= 0.877877;  rate[16][1]= 0.204109;  rate[16][2]= 0.471268;  rate[16][3]= 0.178197;  
    rate[16][4]= 0.079511;  rate[16][5]= 0.248992;  rate[16][6]= 0.321028;  rate[16][7]= 0.136266;  
    rate[16][8]= 0.101117;  rate[16][9]= 0.376588;  rate[16][10]= 0.243227; rate[16][11]= 0.446646; 
    rate[16][12]= 0.184609; rate[16][13]= 0.08587;  rate[16][14]= 0.207143; rate[16][15]= 1.767766; 
    rate[17][0]= 0.030309;  rate[17][1]= 0.046417;  rate[17][2]= 0.010459;  rate[17][3]= 0.011393;  
    rate[17][4]= 0.007732;  rate[17][5]= 0.021248;  rate[17][6]= 0.018844;  rate[17][7]= 0.02399;   
    rate[17][8]= 0.020009;  rate[17][9]= 0.034954;  rate[17][10]= 0.083439; rate[17][11]= 0.023321; 
    rate[17][12]= 0.022019; rate[17][13]= 0.12805;  rate[17][14]= 0.014584; rate[17][15]= 0.035933; 
    rate[17][16]= 0.020437; rate[18][0]= 0.087061;  rate[18][1]= 0.09701;   rate[18][2]= 0.093268;  
    rate[18][3]= 0.051664;  rate[18][4]= 0.042823;  rate[18][5]= 0.062544;  rate[18][6]= 0.0552;    
    rate[18][7]= 0.037568;  rate[18][8]= 0.286027;  rate[18][9]= 0.086237;  rate[18][10]= 0.189842; 
    rate[18][11]= 0.068689; rate[18][12]= 0.073223; rate[18][13]= 0.898663; rate[18][14]= 0.032043; 
    rate[18][15]= 0.121979; rate[18][16]= 0.094617; rate[18][17]= 0.124746; rate[19][0]= 1.230985;  
    rate[19][1]= 0.113146;  rate[19][2]= 0.049824;  rate[19][3]= 0.048769;  rate[19][4]= 0.163831;  
    rate[19][5]= 0.112027;  rate[19][6]= 0.205868;  rate[19][7]= 0.082579;  rate[19][8]= 0.068575;  
    rate[19][9]= 3.65443;   rate[19][10]= 1.337571; rate[19][11]= 0.144587; rate[19][12]= 0.307309; 
    rate[19][13]= 0.247329; rate[19][14]= 0.129315; rate[19][15]= 0.1277;   rate[19][16]= 0.740372; 
    rate[19][17]= 0.022134; rate[19][18]= 0.125733; 
    
       return rate;
    }

    @Override
    public double[] getEmpiricalFrequencies() {
        double[] f = new double[20];
    
    f[0]= 0.078837;         f[1]= 0.051238;         f[2]= 0.042313;         f[3]= 0.053066;         
    f[4]= 0.015175;         f[5]= 0.036713;         f[6]= 0.061924;         f[7]= 0.070852;         
    f[8]= 0.023082;         f[9]= 0.062056;         f[10]= 0.096371;        f[11]= 0.057324;        
    f[12]= 0.023771;        f[13]= 0.043296;        f[14]= 0.043911;        f[15]= 0.063403;        
    f[16]= 0.055897;        f[17]= 0.013272;        f[18]= 0.034399;        f[19]= 0.073101;        
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
