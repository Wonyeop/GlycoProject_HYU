package glycanDB;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

public class WURCS_reader {
	public boolean readWURCS(String filePath) {
		BufferedReader br = null;
		try {
			int br_cnt = 0;
			
			br = new BufferedReader(new FileReader(filePath));
			String line = br.readLine(); // header skip ()
			while((line = br.readLine())!=null) {
				br_cnt++;
				String acc = line.substring(0, line.indexOf(','));
				String wurcs = line.substring(line.indexOf('"')+1, line.length()-1);
				Glycan sugar = new Glycan(acc, wurcs);
				if( br_cnt > 50 ) break;
			}
			br.close();
		} catch (IOException e) {
			System.out.println("I/O exception1");
			e.printStackTrace();
			return false;
		}
		return false;
	}
	
	
	public WURCS_reader() {
		
	}
	
	public static void main(String[] args) {
		WURCS_reader wr = new WURCS_reader();
		wr.readWURCS("C:/Users/user_lee/Desktop/230508_GlyTouCan/glycosmos_glycans_wurcs.csv");
	}
}
