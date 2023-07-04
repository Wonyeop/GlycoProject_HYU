package glycanDB;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

public class Offset_test {
	public static void main(String[] args) throws IOException {
		String line;
		String fileName = "msgf_glycan_3spectrum_test_FDR_converted_modp.mgf";
		
		BufferedReader br = new BufferedReader(new FileReader(fileName));
		long offset = 0;	//	charactors count
		System.out.println("Count offsets");
		int index = 1;
		while( (line=br.readLine()) != null ){
			offset += line.length() +2;
			
			if( line.startsWith("BEGIN IONS") ){
				line = br.readLine();
				offset += line.length() +2;
				line = br.readLine();
				
				while( !line.startsWith("END IONS") ){
					
					String[] spliter = line.split("\\s+");
					if( spliter.length == 2 && Character.isDigit(line.charAt(0)) ){
						try{
							Double.parseDouble(spliter[0]);
							System.out.println(offset);
							break;
						}catch(NumberFormatException e){
							
						}
					}
					
					offset += line.length() +2;
					line = br.readLine();
				}
				offset += line.length() +2;
				
				index++;
			}
			
		}
		br.close();
		
		System.out.println("--------------");
		
		br = new BufferedReader(new FileReader(fileName));
		br.skip(165);
		System.out.println("Scan1");
		System.out.println(br.readLine());
		System.out.println(br.readLine());
		System.out.println(br.readLine());
		br.close();
		br = new BufferedReader(new FileReader(fileName));
		br.skip(14064);
		System.out.println("Scan2");
		System.out.println(br.readLine());
		System.out.println(br.readLine());
		System.out.println(br.readLine());
		br.close();
		br = new BufferedReader(new FileReader(fileName));
		br.skip(17922);
		System.out.println("Scan3");
		System.out.println(br.readLine());
		System.out.println(br.readLine());
		System.out.println(br.readLine());
		br.close();
		
		
		
	}
}
