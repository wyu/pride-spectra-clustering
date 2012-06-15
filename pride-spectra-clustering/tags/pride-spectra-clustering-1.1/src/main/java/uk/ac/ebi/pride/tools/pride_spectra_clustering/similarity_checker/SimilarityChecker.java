package uk.ac.ebi.pride.tools.pride_spectra_clustering.similarity_checker;

import java.util.List;

import uk.ac.ebi.pride.tools.pride_spectra_clustering.util.Peak;

/**
 * Assesses the similarity between two
 * spectra and returns the result of this
 * check as a double. Higher values mean
 * a higher similarity.
 * @author jg
 *
 */
public interface SimilarityChecker {
	/**
	 * Assesses the similarity between the two passed
	 * spectra and returns the result of this assessment
	 * as a double. A higher number reflects a higher
	 * grade of similarity.
	 * @param spectrum1 The first spectrum to compare. The list of Peaks MUST be sorted according to intensity.
	 * @param spectrum2 The second spectrum to compare. The list of Peaks MUST be sorted according to intensity.
	 * @param precursorMz1 Spectrum 1's precursor's m/z
	 * @param precursorMz2 Spectrum 2's precursor's m/z
	 * @param charge1 Spectrum 1's charge state
	 * @param charge2 Spectrum 2's charge state
	 * @return A score indicating the similarity between the two passed spectra.
	 */
	public double assessSimilarity(List<Peak> spectrum1, List<Peak> spectrum2, Double precursorMz1, Double precursorMz2, Double charge1, Double charge2);
}
