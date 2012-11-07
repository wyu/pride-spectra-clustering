package uk.ac.ebi.pride.tools.fast_spectra_clustering;


import java.io.File;
import java.net.URL;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import junit.framework.TestCase;

import org.junit.Before;
import org.junit.Test;

import uk.ac.ebi.pride.tools.jmzreader.model.Spectrum;
import uk.ac.ebi.pride.tools.mgf_parser.MgfFile;
import uk.ac.ebi.pride.tools.mgf_parser.model.Ms2Query;
import uk.ac.ebi.pride.tools.pride_spectra_clustering.SpectraClustering;
import uk.ac.ebi.pride.tools.pride_spectra_clustering.impl.FrankEtAlClustering;
import uk.ac.ebi.pride.tools.pride_spectra_clustering.util.Peak;
import uk.ac.ebi.pride.tools.pride_spectra_clustering.util.SpectraCluster;

public class FrankEtAlClusteringTest extends TestCase {
	private static final SpectraClustering clustering = new FrankEtAlClustering();
	private static File specFile;
	private static MgfFile mgfFile;
	private static List<Spectrum> spectra;

	@Before
	public void setUp() throws Exception {
		URL testFile = getClass().getClassLoader().getResource("spectra_400.0_4.0.mgf");
		assertNotNull(testFile);
		
		if (specFile == null)
			specFile = new File(testFile.toURI());
		
		assertNotNull(specFile);
		
		if (mgfFile == null)
			mgfFile = new MgfFile(specFile);
		
		assertNotNull(mgfFile);
		
		if (spectra == null) {
			spectra = new ArrayList<Spectrum>(mgfFile.getMs2QueryCount());
			Iterator<Ms2Query> it = mgfFile.getMs2QueryIterator();
			while (it.hasNext()) {
				Ms2Query query = it.next();
				if (query.getPrecursorIntensity() == null)
					query.setPeptideIntensity(1.0);
				
				spectra.add(query);
			}
		}
	}
	
	@Test
	public void testClustering1() {
		// set the clustering parameters
		clustering.setClusteringRounds(2);
		clustering.setSimilarityThreshold(0.7);
		
		// do the clustering
		long start = System.currentTimeMillis();
		List<SpectraCluster> generatedCluster = clustering.clusterSpectra(spectra);
		long stop = System.currentTimeMillis();
		
		System.out.println("Clustering done in " + (stop - start) + " msec");
		
		assertEquals(142, generatedCluster.size());
		
		for (int i = 0; i < 3; i++) {
			SpectraCluster cluster = generatedCluster.get(i);
			
			if (i == 0) {
				assertEquals(402.01793764705883, cluster.getAverageMz());
                assertEquals(8, cluster.getClusterSize());
                //System.out.println("Here1");
				boolean peakFound = false;
				for (Peak p : cluster.getConsensusSpectrum()) {
					if (p.getMz() == 764.5533508019715) {
						assertEquals(0.013422155963675658, p.getIntensity());
						peakFound = true;
					}
				}
				assertTrue(peakFound);
			}
		}
	}
	
	@Test
	public void testClustering2() {
		// set the clustering parameters
		clustering.setClusteringRounds(2);
		clustering.setSimilarityThreshold(0.8);
		
		// do the clustering
		long start = System.currentTimeMillis();
		List<SpectraCluster> generatedCluster = clustering.clusterSpectra(spectra);
		long stop = System.currentTimeMillis();
		
		System.out.println("Clustering done in " + (stop - start) + " msec");
		
		assertEquals(168, generatedCluster.size());
		
		for (int i = 0; i < 6; i++) {
			SpectraCluster cluster = generatedCluster.get(i);
			
			if (i == 5) {
				assertEquals(401.67999, cluster.getAverageMz());
				assertEquals(1, cluster.getClusterSize());
				
				boolean peakFound = false;
				for (Peak p : cluster.getConsensusSpectrum()) {
					if (p.getMz() == 233.92139) {
						assertEquals(0.096907934822912, p.getIntensity());
						peakFound = true;
					}
				}
				assertTrue(peakFound);
			}
		}
	}

}
