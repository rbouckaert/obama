package beast.evolution.substitutionmodel;

import beast.evolution.datatype.Aminoacid;
import beast.evolution.datatype.DataType;
/** model data from codonPHYML, which is based on PHYML **/

public class OBAMA_MtArt extends EmpiricalSubstitutionModel {
    @Override
    double[][] getEmpiricalRates() {
        double[][] rate = new double[20][20];

    /* 
     Federico Abascal, April 2005 (c).
     
     This model has been derived from 36 artropoda mitochondrial genomes.
     
     Each gene of the given species was aligned individually. Then, alignments of the whole set 
     of 13 genes where concatenated and passed through GBlocks (Castresana, 2000, in JME) with 
     parameters and output:
     
     Minimum Number Of Sequences For A Conserved Position: 20
     Minimum Number Of Sequences For A Flanking Position: 32
     Maximum Number Of Contiguous Nonconserved Positions: 8
     Minimum Length Of A Block: 10
     Allowed Gap Positions: With Half
     Use Similarity Matrices: Yes
     
     Flank positions of the 40 selected block(s)
     Flanks: [6  22]  [26  44]  [61  70]  [77  143]  [145  185]  [208  236]  [309  640]  
     [644  802]  [831  941]  [956  966]  [973  1062]  [1085  1339]  [1343  1702]  
     [1754  1831]  [1840  1911]  [1916  1987]  [2011  2038]  [2097  2118]  [2125  2143]  
     [2179  2215]  [2243  2268]  [2277  2288]  [2333  2347]  [2476  2518]  [2539  2558]  
     [2600  2613]  [2637  2672]  [2738  2759]  [2784  2839]  [2882  2924]  [2948  3097]  
     [3113  3123]  [3210  3235]  [3239  3322]  [3348  3392]  [3406  3526]  [3588  3617]  
     [3660  3692]  [3803  3830]  [3909  3928]  
     
     New number of positions in MtArt-strict.phy.fasta-gb: <b> 2664 </b> (67% of the original 3933 positions)
     
     The species included in the analysis were:
     Harpiosquilla harpax          [NCBI_TaxID 287944]    
     Ixodes uriae                  [NCBI_TaxID 59655]     
     Heptathela hangzhouensis      [NCBI_TaxID 216259]    
     Triops longicaudatus          [NCBI_TaxID 58777]     
     Gryllotalpa orientalis        [NCBI_TaxID 213494]    
     lepidopsocid RS-2001          [NCBI_TaxID 159971]    
     Locusta migratoria            [NCBI_TaxID 7004]      
     Drosophila yakuba             [NCBI_TaxID 7245]      
     Ostrinia furnacalis           [NCBI_TaxID 93504]     
     Megabalanus volcano           [NCBI_TaxID 266495]    
     Periplaneta fuliginosa        [NCBI_TaxID 36977]     
     Thermobia domestica           [NCBI_TaxID 89055]     
     Aleurochiton aceris           [NCBI_TaxID 266942]    
     Schizaphis graminum           [NCBI_TaxID 13262]     
     Pteronarcys princeps          [NCBI_TaxID 285953]    
     Aleurodicus dugesii           [NCBI_TaxID 30099]     
     Pollicipes polymerus          [NCBI_TaxID 36137]     
     Gomphiocephalus hodgsoni      [NCBI_TaxID 221270]    
     Habronattus oregonensis       [NCBI_TaxID 130930]    
     Speleonectes tulumensis       [NCBI_TaxID 84346]     
     Hutchinsoniella macracantha   [NCBI_TaxID 84335]     
     Haemaphysalis flava           [NCBI_TaxID 181088]    
     Scutigera coleoptrata         [NCBI_TaxID 29022]     
     Vargula hilgendorfii          [NCBI_TaxID 6674]      
     Tricholepidion gertschi       [NCBI_TaxID 89825]     
     Varroa destructor             [NCBI_TaxID 109461]    
     Bombyx mandarina              [NCBI_TaxID 7092]      
     Thyropygus sp.                [NCBI_TaxID 174155]    
     Tribolium castaneum           [NCBI_TaxID 7070]      
     Pagurus longicarpus           [NCBI_TaxID 111067]    
     Limulus polyphemus            [NCBI_TaxID 6850]      
     Tetrodontophora bielanensis   [NCBI_TaxID 48717]     
     Penaeus monodon               [NCBI_TaxID 6687]      
     Daphnia pulex                 [NCBI_TaxID 6669]      
     Apis mellifera                [NCBI_TaxID 7469]      
     Anopheles gambiae             [NCBI_TaxID 7165]      
     
     The topoLOGy used for inferring the model was:
     (((Daph_pulex,Trio_longi),((((((Aleu_aceri,Aleu_duges),Schi_grami),lepi_RS_20),
     ((((Ostr_furna,Bomb_manda),(Dros_yakub,Anop_gambi)),Apis_melli),Trib_casta)),
     ((Gryl_orien,Locu_migra),(Pter_princ,Peri_fulig))),(Tric_gerts,Ther_domes)),
     (Scut_coleo,Thyr_sp),Varg_hilge,Hutc_macra,((((Ixod_uriae,Haem_flava),Varr_destr),
     (Habr_orego,Hept_hangz)),Limu_polyp),(Poll_polym,Mega_volca),(Gomp_hodgs,Tetr_biela),
     ((Pagu_longi,Pena_monod),Harp_harpa),Spel_tulum));
     
     Note this is not the ML topoLOGy but the consensus one (based on morphoLOGical data, 
     phyLOGenetic reconstruction using nuclear genes, etc). Where relationships are
     not clear, a polytomy was introduced (it contains quite a lot of polytomies!).
     
     The model was estimated using (the great and helpful) Ziheng Yang's Paml software package.
     A four-categorized gamma distribution was used to account for heterogeneity (alpha
     was estimated to be 0.47821). Sites with ambiguity data were taken into account.
     
     If you would like the data related to this matrix, please contact fabascal@uvigo.es.
     Federico Abascal (c)2005.
     */
    rate[1][0] = 0.2;     rate[2][0] = 0.2;     rate[2][1] = 0.2;     rate[3][0] = 0.6;
    rate[3][1] = 4.3;     rate[3][2] = 500.2;   rate[4][0] = 253.5;   rate[4][1] = 35.5;
    rate[4][2] = 98.2;    rate[4][3] = 10.6;    rate[5][0] = 0.2;     rate[5][1] = 154.0;
    rate[5][2] = 261.8;   rate[5][3] = 0.2;     rate[5][4] = 0.2;     rate[6][0] = 0.2;
    rate[6][1] = 0.2;     rate[6][2] = 183.0;   rate[6][3] = 861.8;   rate[6][4] = 0.2;
    rate[6][5] = 261.6;   rate[7][0] = 199.8;   rate[7][1] = 0.2;     rate[7][2] = 120.5;
    rate[7][3] = 12.5;    rate[7][4] = 80.5;    rate[7][5] = 2.6;     rate[7][6] = 43.9;
    rate[8][0] = 0.2;     rate[8][1] = 41.3;    rate[8][2] = 179.5;   rate[8][3] = 0.2;
    rate[8][4] = 12.4;    rate[8][5] = 313.5;   rate[8][6] = 15.2;    rate[8][7] = 0.2;
    rate[9][0] = 25.7;    rate[9][1] = 1.8;     rate[9][2] = 21.3;    rate[9][3] = 6.6;
    rate[9][4] = 63.0;    rate[9][5] = 10.5;    rate[9][6] = 6.8;     rate[9][7] = 2.7;
    rate[9][8] = 0.2;     rate[10][0] = 3.7;    rate[10][1] = 1.8;    rate[10][2] = 12.6;
    rate[10][3] = 1.2;    rate[10][4] = 78.7;   rate[10][5] = 16.3;   rate[10][6] = 1.7;
    rate[10][7] = 1.4;    rate[10][8] = 5.5;    rate[10][9] = 514.5;  rate[11][0] = 0.2;
    rate[11][1] = 208.6;  rate[11][2] = 467.3;  rate[11][3] = 1.7;    rate[11][4] = 0.2;
    rate[11][5] = 349.3;  rate[11][6] = 106.3;  rate[11][7] = 0.2;    rate[11][8] = 0.2;
    rate[11][9] = 3.5;    rate[11][10] = 3.8;   rate[12][0] = 120.6;  rate[12][1] = 5.2;
    rate[12][2] = 78.8;   rate[12][3] = 0.2;    rate[12][4] = 312.3;  rate[12][5] = 67.3;
    rate[12][6] = 0.2;    rate[12][7] = 55.7;   rate[12][8] = 0.2;    rate[12][9] = 514.8;
    rate[12][10] = 885.5; rate[12][11] = 105.6; rate[13][0] = 13.1;   rate[13][1] = 4.7;
    rate[13][2] = 19.7;   rate[13][3] = 0.2;    rate[13][4] = 184.1;  rate[13][5] = 0.2;
    rate[13][6] = 0.2;    rate[13][7] = 0.8;    rate[13][8] = 13.8;   rate[13][9] = 117.9;
    rate[13][10] = 262.6; rate[13][11] = 10.7;  rate[13][12] = 321.6; rate[14][0] = 49.3;
    rate[14][1] = 0.2;    rate[14][2] = 16.5;   rate[14][3] = 0.2;    rate[14][4] = 0.2;
    rate[14][5] = 39.3;   rate[14][6] = 7.9;    rate[14][7] = 0.2;    rate[14][8] = 0.8;
    rate[14][9] = 0.2;    rate[14][10] = 12.2;  rate[14][11] = 16.8;  rate[14][12] = 5.3;
    rate[14][13] = 14.6;  rate[15][0] = 673.0;  rate[15][1] = 2.7;    rate[15][2] = 398.4;
    rate[15][3] = 44.4;   rate[15][4] = 664.2;  rate[15][5] = 52.4;   rate[15][6] = 31.5;
    rate[15][7] = 226.0;  rate[15][8] = 10.6;   rate[15][9] = 7.2;    rate[15][10] = 8.2;
    rate[15][11] = 144.2; rate[15][12] = 111.7; rate[15][13] = 36.1;  rate[15][14] = 86.5;
    rate[16][0] = 243.9;  rate[16][1] = 0.2;    rate[16][2] = 165.9;  rate[16][3] = 0.2;
    rate[16][4] = 182.8;  rate[16][5] = 43.7;   rate[16][6] = 43.4;   rate[16][7] = 0.2;
    rate[16][8] = 18.6;   rate[16][9] = 203.7;  rate[16][10] = 47.8;  rate[16][11] = 69.5;
    rate[16][12] = 288.6; rate[16][13] = 13.5;  rate[16][14] = 46.8;  rate[16][15] = 660.4;
    rate[17][0] = 0.2;    rate[17][1] = 0.2;    rate[17][2] = 7.7;    rate[17][3] = 0.2;
    rate[17][4] = 21.6;   rate[17][5] = 6.7;    rate[17][6] = 11.0;   rate[17][7] = 1.9;
    rate[17][8] = 0.2;    rate[17][9] = 0.2;    rate[17][10] = 21.1;  rate[17][11] = 16.0;
    rate[17][12] = 70.7;  rate[17][13] = 53.7;  rate[17][14] = 0.2;   rate[17][15] = 2.4;
    rate[17][16] = 0.2;   rate[18][0] = 1.2;    rate[18][1] = 3.9;    rate[18][2] = 251.2;
    rate[18][3] = 0.2;    rate[18][4] = 72.0;   rate[18][5] = 86.7;   rate[18][6] = 7.7;
    rate[18][7] = 8.6;    rate[18][8] = 191.4;  rate[18][9] = 12.3;   rate[18][10] = 19.8;
    rate[18][11] = 117.1; rate[18][12] = 70.9;  rate[18][13] = 791.6; rate[18][14] = 18.4;
    rate[18][15] = 30.5;  rate[18][16] = 46.0;  rate[18][17] = 37.7;  rate[19][0] = 339.9;
    rate[19][1] = 0.2;    rate[19][2] = 22.6;   rate[19][3] = 0.2;    rate[19][4] = 350.4;
    rate[19][5] = 0.2;    rate[19][6] = 13.6;   rate[19][7] = 2.6;    rate[19][8] = 0.2;
    rate[19][9] = 1854.5; rate[19][10] = 84.7;  rate[19][11] = 26.1;  rate[19][12] = 281.3;
    rate[19][13] = 51.9;  rate[19][14] = 31.7;  rate[19][15] = 60.6;  rate[19][16] = 544.1;
    rate[19][17] = 0.2;   rate[19][18] = 1.6;
    
       for (int i = 0; i < 20; i++) for (int j = i + 1; j < 20; j++) rate[i][j] = rate[j][i]; return rate;
    }

    @Override
    public double[] getEmpiricalFrequencies() {
        double[] f = new double[20];
    
    f[0] = 0.054116;       f[1] = 0.018227;       f[2] = 0.039903;       f[3] = 0.020160;       f[4] = 0.009709;
    f[5] = 0.018781;       f[6] = 0.024289;       f[7] = 0.068183;       f[8] = 0.024518;       f[9] = 0.092639;
    f[10] = 0.148658;      f[11] = 0.021718;      f[12] = 0.061453;      f[13] = 0.088668;      f[14] = 0.041825; // was 0.041826 but does not add to 1
    f[15] = 0.091030;      f[16] = 0.049194;      f[17] = 0.029786;      f[18] = 0.039443;      f[19] = 0.057700; // was 0.057701 but does not add to 1
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
