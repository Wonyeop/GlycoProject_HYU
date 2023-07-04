package datautils;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

public class CountSpectrum {
	public CountSpectrum() {
		
	}
	
	public void process(String file_dir) {
		try {
			File fpath = new File(file_dir);
			if (fpath.isDirectory()) {
				File[] list_input_mgf_files = fpath.listFiles();
				for (int i = 0; i < list_input_mgf_files.length; i++) {
					mgf_reader(list_input_mgf_files[i]);
				}
			} else
				mgf_reader(fpath);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}	
	
	private void mgf_reader(File fmgf) throws IOException{
		BufferedReader br = new BufferedReader(new FileReader(fmgf));
		String line = null;
		
		int cnt_scan = 0;
		while( (line=br.readLine()) != null ) {
			if( line.startsWith("END IONS") ) {
				cnt_scan++;
			}
		}
		System.out.println(fmgf.getName() + "\t" + cnt_scan);
		br.close();
	}
	
	public static void main(String[] args) {
		CountSpectrum cs = new CountSpectrum();
		cs.process("D:\\Glycans\\230313_i-GPA_test\\peak");
	}
}
