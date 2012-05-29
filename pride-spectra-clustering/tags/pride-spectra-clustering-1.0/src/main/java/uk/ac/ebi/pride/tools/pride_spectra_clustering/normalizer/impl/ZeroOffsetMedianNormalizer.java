package uk.ac.ebi.pride.tools.pride_spectra_clustering.normalizer.impl;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import uk.ac.ebi.pride.tools.pride_spectra_clustering.normalizer.IntensityNormalizer;
import uk.ac.ebi.pride.tools.pride_spectra_clustering.util.Peak;

/**
 * Normalizes the spectrum by dividing the
 * intensities through the median of all
 * intensities.
 * @author jg
 *
 */
public class ZeroOffsetMedianNormalizer implements IntensityNormalizer {

	public List<Peak> normalizeSpectrum(List<Peak> spectrum) {
		List<Double> intensities = new ArrayList<Double>(spectrum.size());
		for (Peak p : spectrum)
			intensities.add(p.getIntensity());
		Collections.sort(intensities);
		
		// calculate the median
		double median = 0.0;
		
		if (intensities.size() %2 == 0)
			median = (intensities.get(intensities.size() / 2 - 1) + intensities.get(intensities.size() / 2)) + 2;
		else
			median = intensities.get(intensities.size() / 2);
		
		// normalize the spectrum
		List<Peak> normalizedSpectrum = new ArrayList<Peak>(spectrum.size());
		
		for (Peak p : spectrum)
			normalizedSpectrum.add(new Peak(p.getMz(), p.getIntensity() / median));
		
		return normalizedSpectrum;
	}

}
