package uk.ac.ebi.pride.tools.pride_spectra_clustering.util;

import java.util.Comparator;

public class PeakIntensityComparator  implements Comparator<Peak> {
	private static PeakIntensityComparator instance = new PeakIntensityComparator();
	
	private PeakIntensityComparator() {
		
	}
	
	public static PeakIntensityComparator getInstance() {		
		return instance;
	}

	public int compare(Peak peak1, Peak peak2) {
		if (peak1 == null && peak2 != null)
			return -1;
		if (peak1 == null && peak2 == null)
			return 0;
		if (peak1 != null && peak2 == null)
			return 1;
		
		if (peak1.getIntensity() < peak2.getIntensity())
			return -1;
		if (peak1.getIntensity() > peak2.getIntensity())
			return 1;
		
		return 0;
	}
}
