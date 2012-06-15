package uk.ac.ebi.pride.tools.pride_spectra_clustering.quality_check;

import java.util.List;

import uk.ac.ebi.pride.tools.pride_spectra_clustering.util.Peak;

/**
 * Checks the quality of a given
 * spectrum and returns this quality
 * as a double.
 * @author jg
 *
 */
public interface QualityChecker {
	/**
	 * Assesses the quality of the given spectrum
	 * and returns this assessment as a double. Higher
	 * quality scores should represent a higher
	 * quality of the spectrum.
	 * @param spectrum The spectrum as a map with the m/z values as key and their intensities as value.
	 * @return The spectrum's quality score.
	 */
	public double assessQuality(List<Peak> spectrum);
	/**
	 * Indicates whether the given quality
	 * checker required normalized data.
	 * @return
	 */
	public boolean requiresNormalization();
}
