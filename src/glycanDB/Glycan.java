package glycanDB;

import java.util.Map;

public class Glycan {
	public String accession_id = null;
	private String original_wurcs = null;
	private String original_IUPAC_Condensed = null;
	public String composition = null;
//	public Map<Double> edges = null;
	
	public Glycan(String acc_id, String wurcs) {
		accession_id = acc_id;
		original_wurcs = wurcs;
//		if(!parseWURCS()) System.err.println(-1);
	}
	
	public String getWURCS() {return original_wurcs;}
	public String getGlyTouCanAccession() {return accession_id;}
	
//	private boolean parseDB() {
//	
//	}
	
	public class Node {
		public double mass;
		
	}
}

