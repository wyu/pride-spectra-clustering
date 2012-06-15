package uk.ac.ebi.pride.tools.pride_spectra_clustering.normalizer.impl;

import java.util.ArrayList;
import java.util.List;

import uk.ac.ebi.pride.tools.pride_spectra_clustering.normalizer.IntensityNormalizer;
import uk.ac.ebi.pride.tools.pride_spectra_clustering.util.Peak;

/**
 * Normalizes a spectrum by dividing the intensities
 * through the mean of all intensities in the spectrum.
 * 
 * TODO: add reference
 * @author jg
 *
 */
public class ZeroOffsetMeanNormalizer implements IntensityNormalizer {

	public List<Peak> normalizeSpectrum(List<Peak> spectrum) {
		// calculate the total intensity
		Double totalIntensity = 0.0;
		
		for (Peak p : spectrum)
			totalIntensity += p.getIntensity();
		
		double averageIntensity = totalIntensity / spectrum.size();
		
		// normalize the spectrum
		List<Peak> normalizedSpectrum = new ArrayList<Peak>(spectrum.size());
		
		for (Peak p : spectrum)
			normalizedSpectrum.add(new Peak(p.getMz(), p.getIntensity() / averageIntensity));
		
		return normalizedSpectrum;
	}

}
