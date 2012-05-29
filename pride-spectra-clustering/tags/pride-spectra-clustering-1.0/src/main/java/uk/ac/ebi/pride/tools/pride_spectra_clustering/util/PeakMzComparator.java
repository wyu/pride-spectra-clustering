package uk.ac.ebi.pride.tools.pride_spectra_clustering.util;

import java.util.Comparator;

public class PeakMzComparator implements Comparator<Peak> {
	private static PeakMzComparator instance = new PeakMzComparator();
	
	private PeakMzComparator() {
		
	}
	
	public static PeakMzComparator getInstance() {
		return instance;
	}

	public int compare(Peak peak1, Peak peak2) {
		if (peak1 == null && peak2 != null)
			return -1;
		if (peak1 == null && peak2 == null)
			return 0;
		if (peak1 != null && peak2 == null)
			return 1;
		
		if (peak1.getMz() < peak2.getMz())
			return -1;
		if (peak1.getMz() > peak2.getMz())
			return 1;
		
		return 0;
	}
}
