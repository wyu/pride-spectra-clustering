package uk.ac.ebi.pride.tools.pride_spectra_clustering.normalizer.impl;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import uk.ac.ebi.pride.tools.pride_spectra_clustering.normalizer.IntensityNormalizer;
import uk.ac.ebi.pride.tools.pride_spectra_clustering.util.Peak;

/**
 * Normalizes a spectrum by subtracting
 * the min(intensity) from the current intensity
 * and dividing the result by the range of intensities.
 * 
 * @author jg
 *
 */
public class LinearNormalizer implements IntensityNormalizer {

	public List<Peak> normalizeSpectrum(List<Peak> spectrum) {
		// get the intensities
		List<Double> intensities = new ArrayList<Double>(spectrum.size());
		for (Peak p : spectrum)
			intensities.add(p.getIntensity());
		
		double minIntensity = Collections.min(intensities);
		double maxIntensity = Collections.max(intensities);
		
		// normalize the spectrum
		List<Peak> normalizedSpectrum = new ArrayList<Peak>(spectrum.size());
		
		for (Peak p : spectrum)
			normalizedSpectrum.add(new Peak(p.getMz(), (p.getIntensity() - minIntensity) / (maxIntensity - minIntensity)));
		
		return normalizedSpectrum;
	}

}
