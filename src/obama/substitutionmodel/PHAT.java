package obama.substitutionmodel;

import beast.base.core.Citation;
import beast.base.core.Description;

//from http://www.genome.jp/aaindex/
//	H NGPC000101
//	D Substitution matrix (PHAT) built from hydrophobic and transmembrane regions
//	  of the Blocks database (Ng et al., 2000)
//	R PMID:11108698
//	A Ng, P.C., Henikoff, J.G. and Henikoff, S.
//	T PHAT: a transmembrane-specific substitution matrix
//	J Bioinformatics 16, 760-766 (2000)
//	M rows = ARNDCQEGHILKMFPSTWYV, cols = ARNDCQEGHILKMFPSTWYV


@Description("PHAT: (predicted  hydrophobic and transmembrane) a transmembrane-specific substitution matrix")
@Citation("Ng, Pauline C., Jorja G. Henikoff, and Steven Henikoff. \"PHAT: a transmembrane-specific substitution matrix.\" Bioinformatics 16.9 (2000): 760-766.")
public class PHAT extends ScoreBasedSubstitutionModel {

	@Override
	public int[][] getScores() {
		int [][] scores;
		// from http://www.genome.jp/aaindex/ appears to have double the diagonal
		/* scores = new int[][]{
			{10, -6, -2, -5, 1, -3, -5, 1, -3, 0, -1, -7, -1, -1, -3, 2, 0, -4, -3, 1},
			{-6, 18, -3, -7, -8, -2, -6, -5, -4, -6, -6, -1, -6, -7, -7, -6, -6, -7, -6, -7},
			{-2, -3, 22, 2, -2, 2, 0, -1, 4, -3, -3, -2, -2, -1, -4, 1, -1, -5, 2, -3},
			{-5, -7, 2, 24, -7, 0, 6, -2, -1, -5, -5, -5, -5, -5, -5, -4, -5, -7, -4, -5},
			{1, -8, -2, -7, 14, -5, -7, -2, -7, -3, -2, -10, -2, 0, -8, 1, -1, -4, -1, -2},
			{-3, -2, 2, 0, -5, 18, 1, -2, 2, -3, -3, -1, -1, -2, -3, -1, -3, 1, 0, -3},
			{-5, -6, 0, 6, -7, 1, 24, -3, -1, -5, -5, -4, -5, -5, -5, -3, -5, -7, -2, -5},
			{1, -5, -1, -2, -2, -2, -3, 18, -4, -2, -2, -5, -1, -2, -3, 1, -1, -5, -3, -2},
			{-3, -4, 4, -1, -7, 2, -1, -4, 22, -5, -4, -5, -4, -2, -6, -2, -4, -3, 3, -5},
			{0, -6, -3, -5, -3, -3, -5, -2, -5, 10, 2, -7, 3, 0, -4, -2, -1, -4, -3, 3},
			{-1, -6, -3, -5, -2, -3, -5, -2, -4, 2, 8, -7, 2, 1, -5, -2, -1, -3, -2, 1},
			{-7, -1, -2, -5, -10, -1, -4, -5, -5, -7, -7, 10, -6, -7, -4, -5, -6, -8, -4, -8},
			{-1, -6, -2, -5, -2, -1, -5, -1, -4, 3, 2, -6, 12, 0, -5, -2, 0, -4, -2, 1},
			{-1, -7, -1, -5, 0, -2, -5, -2, -2, 0, 1, -7, 0, 12, -5, -2, -2, 0, 4, -1},
			{-3, -7, -4, -5, -8, -3, -5, -3, -6, -4, -5, -4, -5, -5, 26, -3, -4, -6, -5, -4},
			{2, -6, 1, -4, 1, -1, -3, 1, -2, -2, -2, -5, -2, -2, -3, 12, 1, -5, -2, -2},
			{0, -6, -1, -5, -1, -3, -5, -1, -4, -1, -1, -6, 0, -2, -4, 1, 6, -7, -3, 0},
			{-4, -7, -5, -7, -4, 1, -7, -5, -3, -4, -3, -8, -4, 0, -6, -5, -7, 22, 1, -4},
			{-3, -6, 2, -4, -1, 0, -2, -3, 3, -3, -2, -4, -2, 4, -5, -2, -3, 1, 22, -3},
			{1, -7, -3, -5, -2, -3, -5, -2, -5, 3, 1, -8, 1, -1, -4, -2, 0, -4, -3, 8}};
			*/
		// matrix from the paper = matrix from ttp://www.genome.jp/aaindex/ with diagonal halved
		scores = new int[][]{{5, -6, -2, -5, 1, -3, -5, 1, -3, 0, -1, -7, -1, -1, -3, 2, 0, -4, -3, 1},
				{-6, 9, -3, -7, -8, -2, -6, -5, -4, -6, -6, -1, -6, -7, -7, -6, -6, -7, -6, -7},
				{-2, -3, 11, 2, -2, 2, 0, -1, 4, -3, -3, -2, -2, -1, -4, 1, -1, -5, 2, -3},
				{-5, -7, 2, 12, -7, 0, 6, -2, -1, -5, -5, -5, -5, -5, -5, -4, -5, -7, -4, -5},
				{1, -8, -2, -7, 7, -5, -7, -2, -7, -3, -2, -10, -2, 0, -8, 1, -1, -4, -1, -2},
				{-3, -2, 2, 0, -5, 9, 1, -2, 2, -3, -3, -1, -1, -2, -3, -1, -3, 1, 0, -3},
				{-5, -6, 0, 6, -7, 1, 12, -3, -1, -5, -5, -4, -5, -5, -5, -3, -5, -7, -2, -5},
				{1, -5, -1, -2, -2, -2, -3, 9, -4, -2, -2, -5, -1, -2, -3, 1, -1, -5, -3, -2},
				{-3, -4, 4, -1, -7, 2, -1, -4, 11, -5, -4, -5, -4, -2, -6, -2, -4, -3, 3, -5},
				{0, -6, -3, -5, -3, -3, -5, -2, -5, 5, 2, -7, 3, 0, -4, -2, -1, -4, -3, 3},
				{-1, -6, -3, -5, -2, -3, -5, -2, -4, 2, 4, -7, 2, 1, -5, -2, -1, -3, -2, 1},
				{-7, -1, -2, -5, -10, -1, -4, -5, -5, -7, -7, 5, -6, -7, -4, -5, -6, -8, -4, -8},
				{-1, -6, -2, -5, -2, -1, -5, -1, -4, 3, 2, -6, 6, 0, -5, -2, 0, -4, -2, 1},
				{-1, -7, -1, -5, 0, -2, -5, -2, -2, 0, 1, -7, 0, 6, -5, -2, -2, 0, 4, -1},
				{-3, -7, -4, -5, -8, -3, -5, -3, -6, -4, -5, -4, -5, -5, 13, -3, -4, -6, -5, -4},
				{2, -6, 1, -4, 1, -1, -3, 1, -2, -2, -2, -5, -2, -2, -3, 6, 1, -5, -2, -2},
				{0, -6, -1, -5, -1, -3, -5, -1, -4, -1, -1, -6, 0, -2, -4, 1, 3, -7, -3, 0},
				{-4, -7, -5, -7, -4, 1, -7, -5, -3, -4, -3, -8, -4, 0, -6, -5, -7, 11, 1, -4},
				{-3, -6, 2, -4, -1, 0, -2, -3, 3, -3, -2, -4, -2, 4, -5, -2, -3, 1, 11, -3},
				{1, -7, -3, -5, -2, -3, -5, -2, -5, 3, 1, -8, 1, -1, -4, -2, 0, -4, -3, 4}};
		return scores;
	}

}
