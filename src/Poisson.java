import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;

public class Poisson {
	
	private final double eps=1., dx = 1., dt = 0.1;
	private int N=100;
	private int[] iPlus = new int[N], iMinus = new int[N];
	private double[][][] phi = new double[N][N][N], phiNew = new double[N][N][N];
	private double[][][] rho = new double[N][N][N];
	private double[] E_ijk = new double[3];
	
	public Poisson(int N, double phi0) {
		this.N = N;
		this.init(phi0);
	}
	
	private void setPlusMinus(){
		iPlus = new int[N];
		iMinus = new int[N];
		for(int i=0; i<N; i++){
			iPlus[i]=i+1;
			iMinus[i]=i-1;
		}
		iPlus[N-1]=0;
		iMinus[0]=N-1;
	}
	
	private void init(double phi0){
		
		this.setPlusMinus();
		phi = new double[N][N][N];
		rho = new double[N][N][N];
		
		for(int i=0; i<N; i++)
			for(int j=0; j<N; j++)
				for(int k=0; k<N; k++){
					if (i*j*k < 8) rho[i][j][k] = 1.0;
					else rho[i][j][k] = 0.;
					phi[i][j][k] = phi0;
				}
		
	}
	
	private void getE_ijk(int i, int j, int k){
		E_ijk[0] = -(phi[iPlus[i]][j][k]-phi[iMinus[i]][j][k])/(2.*dx); 
		E_ijk[1] = -(phi[i][iPlus[j]][k]-phi[i][iMinus[j]][k])/(2.*dx); 
		E_ijk[2] = -(phi[i][j][iPlus[k]]-phi[i][j][iMinus[k]])/(2.*dx);
	}
	
	private void setNew(){
		
		phiNew = new double[N][N][N];
		
		for(int i=0; i<N; i++)
			for(int j=0; j<N; j++)
				for(int k=0; k<N; k++)
					phiNew[i][j][k] = (phi[iPlus[i]][j][k]+phi[iMinus[i]][j][k]+
							phi[i][iPlus[j]][k]+phi[i][iMinus[j]][k]+
							phi[i][j][iPlus[k]]+phi[i][j][iMinus[k]])/6.0;
		
		phi = phiNew;
	}
	
	public void evolve(String phiFile) 
			throws FileNotFoundException, UnsupportedEncodingException{
		
		//for(int t=0; t<nt; t++){
			//if(t%dataPoints==0 && t!=0){
				PrintWriter writer1 = new PrintWriter(phiFile, "UTF-8");
				for(int i=0; i<N; i++){
					for(int j=0; j<N; j++)
						for(int k=0; k<N; k++){
							this.getE_ijk(i, j, k);
							writer1.println(i + " " + j + " " + k +  " " + 
									rho[i][j][k] + " " + phi[i][j][k] + " " +
									E_ijk[0] + " " + E_ijk[1] + " " + E_ijk[2]);
						}
					writer1.println();
				}
				writer1.close();
			//}
			//this.setNew();
		//}
		
	}
	
	public static void main(String[] args) 
			throws FileNotFoundException, UnsupportedEncodingException {
		System.out.println("HOLA");
		Poisson poiss = new Poisson(100, 0.0);
		poiss.evolve("Poisson.dat");
		/*Poisson poiss = new Poisson(100, 0.0);
		int dataPoints = 1000;
		String outPhi = "CHphiDensity.dat", 
				outFree = "CHfreeEnergyDensity.dat";
		poiss.evolve(outPhi, outFree, dataPoints, 10000000);*/
	}
	
}
