package test.obama.likelihood;


import org.junit.jupiter.api.Test;

import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.alignment.Sequence;
import beast.base.evolution.branchratemodel.StrictClockModel;
import beast.base.evolution.datatype.UserDataType;
import beast.base.evolution.likelihood.GenericTreeLikelihood;
import beast.base.evolution.substitutionmodel.BinaryCovarion;
import beast.base.evolution.substitutionmodel.Blosum62;
import beast.base.evolution.substitutionmodel.CPREV;
import beast.base.evolution.substitutionmodel.Dayhoff;
import beast.base.evolution.substitutionmodel.Frequencies;
import beast.base.evolution.substitutionmodel.GTR;
import beast.base.evolution.substitutionmodel.GeneralSubstitutionModel;
import beast.base.evolution.substitutionmodel.HKY;
import beast.base.evolution.substitutionmodel.JTT;
import beast.base.evolution.substitutionmodel.JukesCantor;
import beast.base.evolution.substitutionmodel.MTREV;
import beast.base.evolution.substitutionmodel.MutationDeathModel;
import beast.base.evolution.substitutionmodel.SubstitutionModel;
import beast.base.evolution.substitutionmodel.WAG;
import beast.base.evolution.tree.Tree;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import obama.likelihood.MixedTreeLikelihood;
import obama.sitemodel.MixedSiteModel;

import static org.junit.jupiter.api.Assertions.assertEquals;
import test.beast.BEASTTestCase;
import test.beast.evolution.alignment.UncertainAlignmentTest;

/**
 * This test is derived from test.beast.evolution.likelihood.TreeLikelihoodTest
 * It contains all tests with a single component.
 * *
 */
public class MixedTreeLikelihoodSingleComponentTest  {

    public MixedTreeLikelihoodSingleComponentTest() {
        super();
    }

    protected MixedTreeLikelihood newTreeLikelihood() {
    	System.setProperty("java.only", "true");
        return new MixedTreeLikelihood();
    }

    @Test
    public void testJC69Likelihood() throws Exception {
        // Set up JC69 model: uniform freqs, kappa = 1, 0 gamma categories
        Alignment data = BEASTTestCase.getAlignment();
        Tree tree = BEASTTestCase.getTree(data);

        JukesCantor JC = new JukesCantor();
        JC.initAndValidate();

        MixedSiteModel mixedSiteModel = new MixedSiteModel();
        IntegerParameter siteModelIndex = new IntegerParameter();
        siteModelIndex.initByName("value", "0", "dimension", 768);
        mixedSiteModel.initByName("component", JC, "siteModelIndex", siteModelIndex);

        GenericTreeLikelihood likelihood = newTreeLikelihood();
        likelihood.initByName("data", data, "tree", tree, "siteModel", mixedSiteModel);
        double logP = 0;
        logP = likelihood.calculateLogP();
        assertEquals(logP, -1992.2056440317247, BEASTTestCase.PRECISION);

        likelihood.initByName("useAmbiguities", true, "data", data, "tree", tree, "siteModel", mixedSiteModel);
        logP = likelihood.calculateLogP();
        assertEquals(logP, -1992.2056440317247, BEASTTestCase.PRECISION);
    }

    @Test
    public void testJC69LikelihoodWithUncertainCharacters() throws Exception {
    	    	    	
    	Alignment data = UncertainAlignmentTest.getAlignment();
    	Alignment data2 = UncertainAlignmentTest.getUncertainAlignment();
    	double[] logL, logL_uncertain;
    	
    	System.out.println("\nTree A:");
    	Tree tree = UncertainAlignmentTest.getTreeA(data2);    	    	
    	logL = testJC69Likelihood(data,tree);
    	logL_uncertain = testJC69Likelihood(data2,tree);
    	double x1 = -11.853202336328778;
    	double x2 = -12.069603116476458;
    	assertEquals(logL[0], x1, BEASTTestCase.PRECISION);    	
    	assertEquals(logL[1], x1, BEASTTestCase.PRECISION);
    	assertEquals(logL_uncertain[0], x1, BEASTTestCase.PRECISION);    	
    	assertEquals(logL_uncertain[1], x2, BEASTTestCase.PRECISION);    	
    	
    	System.out.println("\nTree B:");
    	tree = UncertainAlignmentTest.getTreeB(data2);
    	logL = testJC69Likelihood(data,tree);
    	logL_uncertain = testJC69Likelihood(data2,tree);
    	double x3 = -12.421114302827698;
    	double x4 = -11.62105662310513;
    	assertEquals(logL[0], x3, BEASTTestCase.PRECISION);    	
    	assertEquals(logL[1], x3, BEASTTestCase.PRECISION);
    	assertEquals(logL_uncertain[0], x3, BEASTTestCase.PRECISION);    	
    	assertEquals(logL_uncertain[1], x4, BEASTTestCase.PRECISION);    	    
    	
    	System.out.println("\nTesting alignment doubling:");
    	Alignment data3 = UncertainAlignmentTest.getUncertainAlignmentDoubled();    	    	
    	logL_uncertain = testJC69Likelihood(data3,tree);
    	assertEquals(logL_uncertain[0], 2 * x3, BEASTTestCase.PRECISION);    	
    	assertEquals(logL_uncertain[1], 2 * x4, BEASTTestCase.PRECISION);    	    
    	
    }        
    
    public double[] testJC69Likelihood(Alignment data, Tree tree) throws Exception {
        // Set up JC69 model: uniform freqs, kappa = 1, 0 gamma categories                              
        JukesCantor JC = new JukesCantor();
        JC.initAndValidate();

        MixedSiteModel mixedSiteModel = new MixedSiteModel();
        IntegerParameter siteModelIndex = new IntegerParameter();
        siteModelIndex.initByName("value", "0", "dimension", data.getSiteCount());
        mixedSiteModel.initByName("component", JC, "siteModelIndex", siteModelIndex);

        // NB The rate in the JC model used here is actually alpha * 3 in the usual sense, because
        // it's divided by 3 before multiplying in the exponent (not sure why)

        System.out.println("Without tip likelihoods:");
        GenericTreeLikelihood likelihood = newTreeLikelihood();
        StrictClockModel clockModel = new StrictClockModel();
        clockModel.initByName("clock.rate", "0.6"); 
        likelihood.initByName("data", data, "tree", tree, "siteModel", mixedSiteModel, "scaling", MixedTreeLikelihood.Scaling.none, 
        		"branchRateModel", clockModel);        
        double[] logP = new double[2];
        logP[0] = likelihood.calculateLogP();
        System.out.println(logP[0]);

        System.out.println("With tip likelihoods:");
        likelihood.initByName("useTipLikelihoods", true, "data", data, "tree", tree, "siteModel", mixedSiteModel, "scaling", MixedTreeLikelihood.Scaling.none);
        logP[1]= likelihood.calculateLogP();
        System.out.println(logP[1]);

        return logP;
    }
    
    @Test
    public void testAscertainedJC69Likelihood() throws Exception {
        // as testJC69Likelihood but with ascertained alignment
        Alignment data = BEASTTestCase.getAscertainedAlignment();
        Tree tree = BEASTTestCase.getTree(data);

        Frequencies freqs = new Frequencies();
        freqs.initByName("data", data,
                "estimate", false);

        HKY hky = new HKY();
        hky.initByName("kappa", "1.0", "frequencies", freqs);

        MixedSiteModel mixedSiteModel = new MixedSiteModel();
        IntegerParameter siteModelIndex = new IntegerParameter();
        siteModelIndex.initByName("value", "0", "dimension", data.getSiteCount());
        mixedSiteModel.initByName("component", hky, "siteModelIndex", siteModelIndex);

        GenericTreeLikelihood likelihood = newTreeLikelihood();
        likelihood.initByName("data", data, "tree", tree, "siteModel", mixedSiteModel);

        double logP = 0;
        logP = likelihood.calculateLogP();
        // the following number comes from Beast 1.6        
        assertEquals(logP, -737.7140695360017, BEASTTestCase.PRECISION);
    }

    @Test
    public void testK80Likelihood() throws Exception {
        // Set up K80 model: uniform freqs, kappa = 27.402591, 0 gamma categories
        Alignment data = BEASTTestCase.getAlignment();
        Tree tree = BEASTTestCase.getTree(data);

        Frequencies freqs = new Frequencies();
        freqs.initByName("data", data,
                "estimate", false);

        HKY hky = new HKY();
        hky.initByName("kappa", "27.40259", "frequencies", freqs);

        MixedSiteModel mixedSiteModel = new MixedSiteModel();
        IntegerParameter siteModelIndex = new IntegerParameter();
        siteModelIndex.initByName("value", "0", "dimension", data.getSiteCount());
        mixedSiteModel.initByName("component", hky, "siteModelIndex", siteModelIndex);

        MixedTreeLikelihood likelihood = newTreeLikelihood();
        likelihood.initByName("data", data, "tree", tree, "siteModel", mixedSiteModel);

        double logP = 0;
        logP = likelihood.calculateLogP();
        assertEquals(logP, -1856.303048876734, BEASTTestCase.PRECISION);

        likelihood.initByName("useAmbiguities", true, "data", data, "tree", tree, "siteModel", mixedSiteModel);
        logP = likelihood.calculateLogP();
        assertEquals(logP, -1856.303048876734, BEASTTestCase.PRECISION);
    }

    @Test
    public void testHKY85Likelihood() throws Exception {
        // Set up HKY85 model: estimated freqs, kappa = 29.739445, 0 gamma categories
        Alignment data = BEASTTestCase.getAlignment();
        Tree tree = BEASTTestCase.getTree(data);

        Frequencies freqs = new Frequencies();
        freqs.initByName("data", data);

        HKY hky = new HKY();
        hky.initByName("kappa", "29.739445", "frequencies", freqs);

        MixedSiteModel mixedSiteModel = new MixedSiteModel();
        IntegerParameter siteModelIndex = new IntegerParameter();
        siteModelIndex.initByName("value", "0", "dimension", data.getSiteCount());
        mixedSiteModel.initByName("component", hky, "siteModelIndex", siteModelIndex);

        MixedTreeLikelihood likelihood = newTreeLikelihood();
        likelihood.initByName("data", data, "tree", tree, "siteModel", mixedSiteModel);

        double logP = 0;
        logP = likelihood.calculateLogP();
        assertEquals(logP, -1825.2131708068507, BEASTTestCase.PRECISION);

        likelihood.initByName("useAmbiguities", true, "data", data, "tree", tree, "siteModel", mixedSiteModel);
        logP = likelihood.calculateLogP();
        assertEquals(logP, -1825.2131708068507, BEASTTestCase.PRECISION);
    }


    @Test
    public void testGTRLikelihood() throws Exception {
        // Set up GTR model: no gamma categories, no proportion invariant
        Alignment data = BEASTTestCase.getAlignment();
        Tree tree = BEASTTestCase.getTree(data);

        Frequencies freqs = new Frequencies();
        freqs.initByName("data", data);

        GTR gtr = new GTR();
        gtr.initByName("frequencies", freqs);

        MixedSiteModel mixedSiteModel = new MixedSiteModel();
        IntegerParameter siteModelIndex = new IntegerParameter();
        siteModelIndex.initByName("value", "0", "dimension", data.getSiteCount());
        mixedSiteModel.initByName("component", gtr, "siteModelIndex", siteModelIndex);

        MixedTreeLikelihood likelihood = newTreeLikelihood();
        likelihood.initByName("data", data, "tree", tree, "siteModel", mixedSiteModel);

        double logP = 0;
        logP = likelihood.calculateLogP();
        assertEquals(logP, -1969.145839307625, BEASTTestCase.PRECISION);

        likelihood.initByName("useAmbiguities", false, "data", data, "tree", tree, "siteModel", mixedSiteModel);
        logP = likelihood.calculateLogP();
        assertEquals(logP, -1969.145839307625, BEASTTestCase.PRECISION);
    }

    void aminoacidModelTest(SubstitutionModel substModel, double expectedValue) throws Exception {
        Alignment data = BEASTTestCase.getAminoAcidAlignment();
        Tree tree = BEASTTestCase.getAminoAcidTree(data);
        MixedSiteModel mixedSiteModel = new MixedSiteModel();
        IntegerParameter siteModelIndex = new IntegerParameter();
        siteModelIndex.initByName("value", "0", "dimension", data.getSiteCount());
        mixedSiteModel.initByName("component", substModel, "siteModelIndex", siteModelIndex);

        MixedTreeLikelihood likelihood = newTreeLikelihood();
        likelihood.initByName("data", data, "tree", tree, "siteModel", mixedSiteModel);
        double logP = 0;
        logP = likelihood.calculateLogP();
        assertEquals(expectedValue, logP, BEASTTestCase.PRECISION);
    }

    @Test
    public void testAminoAcidLikelihoodWAG() throws Exception {
        // Set up WAG model
        WAG wag = new WAG();
        wag.initAndValidate();
        aminoacidModelTest(wag, -338.6388785157248);

    }

    @Test
    public void testAminoAcidLikelihoodJTT() throws Exception {
        // JTT
        JTT jtt = new JTT();
        jtt.initAndValidate();
        aminoacidModelTest(jtt, -338.80761792179726);

    }

    @Test
    public void testAminoAcidLikelihoodBlosum62() throws Exception {
        // Blosum62
        Blosum62 blosum62 = new Blosum62();
        blosum62.initAndValidate();
        aminoacidModelTest(blosum62, -345.3825963600176);

    }

    @Test
    public void testAminoAcidLikelihoodDayhoff() throws Exception {
        // Dayhoff
        Dayhoff dayhoff = new Dayhoff();
        dayhoff.initAndValidate();
        aminoacidModelTest(dayhoff, -340.6149187667345);
    }

    @Test
    public void testAminoAcidLikelihoodcpRev() throws Exception {
        // cpRev
        CPREV cpRev = new CPREV();
        cpRev.initAndValidate();
        aminoacidModelTest(cpRev, -348.71458467304154);
    }

    @Test
    public void testAminoAcidLikelihoodMTRev() throws Exception {
        // MTRev
        MTREV mtRev = new MTREV();
        mtRev.initAndValidate();
        aminoacidModelTest(mtRev, -369.4791633617842);

    }


    @Test
    public void testSDolloLikelihood() throws Exception {
        UserDataType dataType = new UserDataType();
        dataType.initByName("states", 2, "codeMap", "0=1, 1=0, ?=0 1, -=0 1");
        Alignment data = new Alignment();

        Sequence German_ST = new Sequence("German_ST", BEASTTestCase.German_ST.dataInput.get());
        Sequence Dutch_List = new Sequence("Dutch_List", BEASTTestCase.Dutch_List.dataInput.get());
        ;
        Sequence English_ST = new Sequence("English_ST", BEASTTestCase.English_ST.dataInput.get());
        ;
        Sequence French = new Sequence("French", BEASTTestCase.French.dataInput.get());
        ;
        Sequence Italian = new Sequence("Italian", BEASTTestCase.Italian.dataInput.get());
        ;
        Sequence Spanish = new Sequence("Spanish", BEASTTestCase.Spanish.dataInput.get());
        ;


        data.initByName("sequence", German_ST, "sequence", Dutch_List, "sequence", English_ST, "sequence", French, "sequence", Italian, "sequence", Spanish,
                "userDataType", dataType
        );

        Tree tree = BEASTTestCase.getTree(data, "((English_ST:0.22743347188019544,(German_ST:0.10557648379843088,Dutch_List:0.10557648379843088):0.12185698808176457):1.5793160946109988,(Spanish:0.11078392189606047,(Italian:0.10119772534558173,French:0.10119772534558173):0.009586196550478737):1.6959656445951337)");

        RealParameter frequencies = new RealParameter("1 0");
        Frequencies freqs = new Frequencies();
        freqs.initByName("frequencies", frequencies);

        RealParameter deathprob = new RealParameter("1.7");
        MutationDeathModel SDollo = new MutationDeathModel();
        SDollo.initByName("deathprob", deathprob, "frequencies", freqs);

        MixedSiteModel mixedSiteModel = new MixedSiteModel();
        IntegerParameter siteModelIndex = new IntegerParameter();
        siteModelIndex.initByName("value", "0", "dimension", data.getSiteCount());
        mixedSiteModel.initByName("component", SDollo, "siteModelIndex", siteModelIndex);

        MixedTreeLikelihood likelihood = newTreeLikelihood();
        likelihood.initByName("data", data, "tree", tree, "siteModel", mixedSiteModel);

        double logP = 0;
        likelihood.initByName("data", data, "tree", tree, "siteModel", mixedSiteModel, "useAmbiguities", true);
        logP = likelihood.calculateLogP();
        // beast1 xml gives -3551.6436
        assertEquals(logP, -3551.6436270344648, BEASTTestCase.PRECISION);
    }


    @Test
    public void testBinaryCovarionLikelihood() throws Exception {
        Alignment data = BEASTTestCase.getCovarionAlignment();
        Tree tree = BEASTTestCase.getTree(data, "((English_ST:0.22743347188019544,(German_ST:0.10557648379843088,Dutch_List:0.10557648379843088):0.12185698808176457):1.5793160946109988,(Spanish:0.11078392189606047,(Italian:0.10119772534558173,French:0.10119772534558173):0.009586196550478737):1.6959656445951337)");


        RealParameter alpha = new RealParameter("0.284");
        RealParameter switchRate = new RealParameter("0.829");
        RealParameter frequencies = new RealParameter("0.683 0.317");
        RealParameter hfrequencies = new RealParameter("0.5 0.5");
        BinaryCovarion covarion = new BinaryCovarion();
        covarion.initByName("alpha", alpha, "switchRate", switchRate, "vfrequencies", frequencies, "hfrequencies", hfrequencies);

        MixedSiteModel mixedSiteModel = new MixedSiteModel();
        IntegerParameter siteModelIndex = new IntegerParameter();
        siteModelIndex.initByName("value", "0", "dimension", data.getSiteCount());
        mixedSiteModel.initByName("component", covarion, "siteModelIndex", siteModelIndex);

        MixedTreeLikelihood likelihood = newTreeLikelihood();
        likelihood.initByName("data", data, "tree", tree, "siteModel", mixedSiteModel, "useSitesNotPatterns", true);

        double logP = 0;
        likelihood.initByName("useAmbiguities", true, "data", data, "tree", tree, "siteModel", mixedSiteModel);
        logP = likelihood.calculateLogP();
        // beast1 xml gives -1730.5363
        assertEquals(logP, -1730.53631739, BEASTTestCase.PRECISION);
    }

    @Test
    public void testMarginalisationOfLikelihoodBinary() throws Exception {
    	// test summation over all patterns adds to 1 for binary data

    	Sequence German_ST =  new Sequence("German_ST", "           10110010");
    	Sequence Dutch_List = new Sequence("Dutch_List", "          11010100");
    	Sequence English_ST = new Sequence("English_ST", "          11101000");

        Alignment data = new Alignment();
        data.initByName("sequence", German_ST, "sequence", Dutch_List, "sequence", English_ST,
                "dataType", "binary"
        );

        Tree tree = BEASTTestCase.getTree(data, "(English_ST:0.22743347188019544,(German_ST:0.10557648379843088,Dutch_List:0.10557648379843088):0.12185698808176457):0.0;");


        RealParameter frequencies = new RealParameter("0.683 0.317");
        Frequencies freqs = new Frequencies();
        freqs.initByName("frequencies", frequencies);
        GeneralSubstitutionModel covarion = new GeneralSubstitutionModel();
        covarion.initByName("frequencies", freqs, "rates", new RealParameter("1.0 1.0"));

        MixedSiteModel mixedSiteModel = new MixedSiteModel();
        IntegerParameter siteModelIndex = new IntegerParameter();
        siteModelIndex.initByName("value", "0", "dimension", data.getSiteCount());
        mixedSiteModel.initByName("component", covarion, "siteModelIndex", siteModelIndex);

        MixedTreeLikelihood likelihood = newTreeLikelihood();
        likelihood.initByName("data", data, "tree", tree, "siteModel", mixedSiteModel);
        
        likelihood.initByName("useAmbiguities", false, "data", data, "tree", tree, "siteModel", mixedSiteModel);
        likelihood.calculateLogP();
        double [] logPs = likelihood.getPatternLogLikelihoods();
        double P = 0;
        for (double d : logPs) {
        	P += Math.exp(d);
        }
        assertEquals(P, 1.0, BEASTTestCase.PRECISION);
    }

    @Test
    public void testMarginalisationOfLikelihoodNucleotide() throws Exception {
    	// test summation over all patterns adds to 1 for nucleotide data

    	Sequence German_ST =  new Sequence("German_ST", "           AAAAAAAAAAAAAAAA CCCCCCCCCCCCCCCC GGGGGGGGGGGGGGGG TTTTTTTTTTTTTTTT");
    	Sequence Dutch_List = new Sequence("Dutch_List", "          AAAACCCCGGGGTTTT AAAACCCCGGGGTTTT AAAACCCCGGGGTTTT AAAACCCCGGGGTTTT");
    	Sequence English_ST = new Sequence("English_ST", "          ACGTACGTACGTACGT ACGTACGTACGTACGT ACGTACGTACGTACGT ACGTACGTACGTACGT");

        Alignment data = new Alignment();
        data.initByName("sequence", German_ST, "sequence", Dutch_List, "sequence", English_ST,
                "dataType", "nucleotide"
        );

        Tree tree = BEASTTestCase.getTree(data, "(English_ST:0.22743347188019544,(German_ST:0.10557648379843088,Dutch_List:0.10557648379843088):0.12185698808176457):0.0;");


        RealParameter frequencies = new RealParameter("0.2 0.3 0.4 0.1");
        Frequencies freqs = new Frequencies();
        freqs.initByName("frequencies", frequencies);
        GeneralSubstitutionModel gsm = new GeneralSubstitutionModel();
        gsm.initByName("rates", "1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0", "frequencies", freqs);

        MixedSiteModel mixedSiteModel = new MixedSiteModel();
        IntegerParameter siteModelIndex = new IntegerParameter();
        siteModelIndex.initByName("value", "0", "dimension", data.getSiteCount());
        mixedSiteModel.initByName("component", gsm, "siteModelIndex", siteModelIndex);

        MixedTreeLikelihood likelihood = newTreeLikelihood();
        likelihood.initByName("data", data, "tree", tree, "siteModel", mixedSiteModel);

        likelihood.initByName("useAmbiguities", false, "data", data, "tree", tree, "siteModel", mixedSiteModel);
        likelihood.calculateLogP();
        double [] logPs = likelihood.getPatternLogLikelihoods();
        double P = 0;
        for (double d : logPs) {
        	P += Math.exp(d);
        }
        assertEquals(P, 1.0, BEASTTestCase.PRECISION);
    }

} // class TreeLikelihoodTest
