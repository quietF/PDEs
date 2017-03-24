import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;

public class Poisson {
	
	private final double eps=1., dx = 1., dt = 0.1;
	private int N=100;
	private int[] iPlus = new int[N+1], iMinus = new int[N+1];
	private double[][][] phi = new double[N+1][N+1][N+1], 
			phiNew = new double[N+1][N+1][N+1];
	private double[][][] rho = new double[N+1][N+1][N+1];
	private double[] E_ijk = new double[3];
	
	public Poisson(int N, double phi0) {
		this.N = N;
		this.init();
	}
	
	private void setPlusMinus(){
		iPlus = new int[N+1];
		iMinus = new int[N+1];
		for(int i=1; i<N; i++){
			iPlus[i]=i+1;
			iMinus[i]=i-1;
		}
	}
	
	private void init(){
		
		this.setPlusMinus();
		phi = new double[N+1][N+1][N+1];
		rho = new double[N+1][N+1][N+1];
		
		for(int i=1; i<N; i++)
			for(int j=1; j<N; j++)
				for(int k=1; k<N; k++){
					if ((i-N/2)*(i-N/2)+(j-N/2)*(j-N/2)+
							(k-N/2)*(k-N/2) < N) rho[i][j][k] = 1.0;
					else rho[i][j][k] = 0.;
				}
		
		this.setNew();
		
	}
	
	private void getE_ijk(int i, int j, int k){
		E_ijk[0] = -(phi[iPlus[i]][j][k]-phi[iMinus[i]][j][k])/(2.*dx); 
		E_ijk[1] = -(phi[i][iPlus[j]][k]-phi[i][iMinus[j]][k])/(2.*dx); 
		E_ijk[2] = -(phi[i][j][iPlus[k]]-phi[i][j][iMinus[k]])/(2.*dx);
	}
	
	private void setNew(){
		
		phiNew = new double[N+1][N+1][N+1];
		
		for(int i=1; i<N; i++)
			for(int j=1; j<N; j++)
				for(int k=1; k<N; k++)
					phiNew[i][j][k] = (phi[iPlus[i]][j][k]+phi[iMinus[i]][j][k]+
							phi[i][iPlus[j]][k]+phi[i][iMinus[j]][k]+
							phi[i][j][iPlus[k]]+phi[i][j][iMinus[k]]+
							rho[i][j][k]/eps)/6.0;
		
		phi = phiNew;
	}
	
	public void evolve(String phiFile) 
			throws FileNotFoundException, UnsupportedEncodingException{
		
		for(int t=0; t<10000000; t++){
			if(t%10000==0 && t!=0){
				PrintWriter writer1 = new PrintWriter(phiFile, "UTF-8");
				for(int i=1; i<N; i++){
					for(int j=1; j<N; j++)
						for(int k=1; k<N; k++){
							this.getE_ijk(i, j, k);
							writer1.println(i + " " + j + " " + k +  " " + 
									rho[i][j][k] + " " + phi[i][j][k] + " " +
									E_ijk[0] + " " + E_ijk[1] + " " + E_ijk[2]);
						}
					writer1.println();
				}
				writer1.close();
			}
			this.setNew();
		}
		
	}
	
	public static void main(String[] args) 
			throws FileNotFoundException, UnsupportedEncodingException {
		System.out.println("HOLA");
		Poisson poiss = new Poisson(50, 0.0);
		poiss.evolve("Poisson.dat");
	}
	
}
