package uk.ac.ebi.pride.tools.pride_spectra_clustering.normalizer.impl;

import java.util.ArrayList;
import java.util.List;

import uk.ac.ebi.pride.tools.pride_spectra_clustering.normalizer.IntensityNormalizer;
import uk.ac.ebi.pride.tools.pride_spectra_clustering.util.Peak;

/**
 * Normalizes a spectrum by subtracting
 * the mean of all intensities from the
 * specific intensity and dividing it
 * through the standard deviation. 
 * 
 * This method was discussed by Borgaonkar et al. (2010)
 * in a comparison of normalization methods for MALDI/SELDI-TOF
 * @author jg
 *
 */
public class ZScoreNormalizer implements IntensityNormalizer {

	public List<Peak> normalizeSpectrum(List<Peak> spectrum) {
		// calculate the total intensity
		Double totalIntensity = 0.0;
		
		for (Peak p : spectrum)
			totalIntensity += p.getIntensity();
		
		double averageIntensity = totalIntensity / spectrum.size();
		
		// calculate the standard deviation
		Double sqrDifferences = 0.0;
		
		for (Peak p : spectrum)
			sqrDifferences += Math.pow(p.getIntensity() - averageIntensity, 2);
		
		Double standardDeviation = Math.sqrt(sqrDifferences / spectrum.size());
		
		// normalize the spectrum
		List<Peak> normalizedSpectrum = new ArrayList<Peak>(spectrum.size());
		
		for (Peak p : spectrum)
			normalizedSpectrum.add(new Peak(p.getMz(), (p.getIntensity() - averageIntensity) / standardDeviation));
		
		return normalizedSpectrum;
	}

}
