package glycanDB;

public class Glycan {
	private String accession_id = null;
	private String original_wurcs = null;
	public double wurcs_ver;
	public int BMU;
	public int MLU;
	public int setM;
	public String composition = null;
	public String structure = null;
	
	public Glycan(String acc_id, String wurcs) {
		accession_id = acc_id;
		original_wurcs = wurcs;
		if(!parseWURCS()) System.err.println(-1);
	}
	
	public String getWURCS() {return original_wurcs;}
	public String getGlyTouCanAccession() {return accession_id;}
	
	private boolean parseWURCS() {
		String subWURCS = original_wurcs;
		System.out.println(original_wurcs);
		int idx_sugar = -1;
		boolean ok_parse = false;
		
		while( (idx_sugar = subWURCS.indexOf('[')) >= 0 ) {
			int idx_end = subWURCS.indexOf(']');
			String monosaccharide = subWURCS.substring(idx_sugar+1,idx_end);
			subWURCS = subWURCS.substring(idx_end+1);
			System.out.println(monosaccharide);
		}
		
		return ok_parse;
	}
}
