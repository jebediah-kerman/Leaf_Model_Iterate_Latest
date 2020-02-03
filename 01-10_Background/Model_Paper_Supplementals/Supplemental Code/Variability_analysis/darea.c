//Source code for temporal variability of the growth in cell area

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main (int argc, char * argv[]){
//Restore the data to arrays
//Col 6:Wild type flower number 6
	type=6;
//time:Time frames N:Number of mesh points CellNout:The label for the outside of the sepal
	int type,time[10],N[10],CellNout[10];
	time[0]=0;time[1]=12;time[2]=24;time[3]=36;time[4]=48;
	N[0]=224081;N[1]=192757;N[2]=275129;N[3]=320182;N[4]=414384;
//loop of time frames
	int t;
	for (t=1; t<=3; t++) { 
	//Read the data files, Create the output file
		int i,j,k,l;
		FILE *file1,*file11,*file2,*file22,*file3,*file4,*file5,*file6,*file7,*file8,*file9,*file10,*file12,*file111,*file222,*file333;
		char filename1[100],filename11[100],filename2[100],filename22[100],filename3[100],filename4[100],filename5[100],filename6[100],filename7[100],filename8[100],filename9[100],filename10[100],filename12[100],filename111[100],filename222[100],filename333[100];
		sprintf(filename1,"%d_%d_mesh.Int.rot.txt",type,time[t]);
		sprintf(filename2,"%d_%d_mesh.near.txt",type,time[t]);
		sprintf(filename11,"%d_%d-%d.dat",type,time[t-1],time[t]);
		sprintf(filename22,"%d_%d-%d.dat",type,time[t],time[t+1]);
		sprintf(filename111,"%d_%d_area.dat",type,time[t-1]);
		sprintf(filename222,"%d_%d_area.dat",type,time[t]);
		sprintf(filename333,"%d_%d_area.dat",type,time[t+1]);
		sprintf(filename4,"%d_%d_darea.dat",type,time[t]);
		file1 = fopen(filename1,"r");		//reading the mesh intensity file
		file2 = fopen(filename2,"r");		//reading the mesh neighbor file
		file11 = fopen(filename11,"r");		//reading lineage file previous-current
		file22 = fopen(filename22,"r");		//reading lineage file current-next
		file111 = fopen(filename111,"r");	//reading area at the previous time frame
		file222 = fopen(filename222,"r");	//reading area at the current time frame
		file333 = fopen(filename333,"r");	//reading area at the next time frame
		file4 = fopen(filename4,"w");		//writing the temporal variability of areal growh file for each mesh
		//reading the mesh intensity file,reading the mesh neighbor file
		double *x,*y,*z;
		int sute,*MTnum,*Cnum,**MTnear;
		x = (double*)calloc(N[t]+10, sizeof(double));
		y = (double*)calloc(N[t]+10, sizeof(double));
		z = (double*)calloc(N[t]+10, sizeof(double));
		Cnum = (int*)calloc(N[t]+10, sizeof(double));
		MTnear = (int**)calloc(20, sizeof(double));
		for (i=1; i<=16; i++)MTnear[i] = (int*)calloc(N[t]+10, sizeof(double));
		for (i=0; i<N[t]; i++) {
			fscanf(file1,"%d%lf%lf%lf%d\n",&sute,&x[i],&y[i],&z[i],&Cnum[i]);
			fscanf(file2,"%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d\n",&MTnum[i],&sute2,&MTnear[1][i],&MTnear[2][i],&MTnear[3][i],&MTnear[4][i],&MTnear[5][i],&MTnear[6][i],&MTnear[7][i],&MTnear[8][i],&MTnear[9][i],&MTnear[10][i],&MTnear[11][i],&MTnear[12][i],&MTnear[13][i],&MTnear[14][i],&MTnear[15][i],&MTnear[16][i]);
		}
		//reading lineage file from previous to current
		int *sun,*parent;
		sun = (int*)calloc(1000, sizeof(double));
		parent = (int*)calloc(1000, sizeof(double));
		for (i=0; i<700; i++) {
			fscanf(file22,"%d%d\n",&sun[i],&parent[i]);
		}
		//reading lineage file from current to next
		int *Sun,*Parent;
		Sun = (int*)calloc(1000, sizeof(double));
		Parent = (int*)calloc(1000, sizeof(double));
		for (i=0; i<700; i++) {
			fscanf(file22,"%d%d\n",&Sun[i],&Parent[i]);
		}
		//reading area at the previous time frame
		int *Cnumareaprev;
		double *areaprev;
		Cnumareaprev = (int*)calloc(10000+10, sizeof(double));
		areaprev = (double*)calloc(10000+10, sizeof(double));
		for (i=0; i<700; i++) {
			fscanf(file111,"%d%lf\n",&Cnumareaprev[i],&areaprev[i]);
		}
		//reading area at the current time frame
		int *Cnumareacur;
		double *areacur;
		Cnumareacur = (int*)calloc(10000+10, sizeof(double));
		areacur = (double*)calloc(10000+10, sizeof(double));
		for (i=0; i<700; i++) {
			fscanf(file222,"%d%lf\n",&Cnumareacur[i],&areacur[i]);
		}
		//reading area at the next time frame
		int *Cnumareanext;
		double *areanext;
		Cnumareanext = (int*)calloc(10000+10, sizeof(double));
		areanext = (double*)calloc(10000+10, sizeof(double));
		for (i=0; i<700; i++) {
			fscanf(file333,"%d%lf\n",&Cnumareanext[i],&areanext[i]);
		}
		fclose(file1);
		fclose(file2);
		fclose(file11);
		fclose(file22);
		fclose(file111);
		fclose(file222);
		fclose(file333);
		//Extract multiple step lineage from previous to next
		int *Sun3,*Parent3;
		Sun3 = (int*)calloc(1000, sizeof(double));
		Parent3 = (int*)calloc(1000, sizeof(double));
		for (j=0; j<700; j++) {
			for (i=0; i<10000; i++) {
				if(i==Sun[j]){
					Sun3[j]=Sun[j];
					for (k=0; k<700; k++) {
						if(Parent[j]==sun[k]){
							Parent3[j]=parent[k];
						}
					}
				}
			}
		}
		//rearrange the order of two area files (prev, curr, next) which are always consistent with previous frame
		int *Cnumareaprev1,*Cnumareacur1,*Cnumareanext1;
		double *areaprev1,*areacur1,*areanext1;
		Cnumareaprev1 = (int*)calloc(10000+10, sizeof(double));
		Cnumareacur1 = (int*)calloc(10000+10, sizeof(double));
		Cnumareanext1 = (int*)calloc(10000+10, sizeof(double));
		areaprev1 = (double*)calloc(10000+10, sizeof(double));
		areacur1 = (double*)calloc(10000+10, sizeof(double));
		areanext1 = (double*)calloc(10000+10, sizeof(double));
		double sumaa,sumbb;
		for (i=0; i<700; i++) {
			sumaa=0.0;
			for (j=0; j<700; j++) {
				//rearrange the orger of current time frame 
				if(sun[j]==Cnumareacur[i] && Cnumareacur[i]>0){
					sumaa=sumaa+areacur[i];
					Cnumareacur1[i]=sun[j];
				}
				for (k=0; k<700; k++) {
					//rearrange the orger of previous time frame 
					if(parent[j]==Cnumareaprev[k] && Cnumareaprev[k]>0){
						areaprev1[i]=areaprev[k];
						Cnumareaprev1[i]=parent[j];
					}
				}
			}
			areacur1[i]=sumaa;
			//rearrange the orger of next time files 
			sumbb=0.0;
			for (j=0; j<700; j++) {
				if(Sun3[j]==Cnumareanext[i] && Cnumareanext[i]>0){
					sumbb=sumbb+areanext[i];
					Cnumareanext1[i]=Sun3[j];
				}
			}
			areanext1[i]=sumbb;
		}
//Calculate temporal variability in the growth of cell area
		double *darea,areag0,areag1;
		double *p_darea,thre_darea,bin1_darea,bin_darea=0.3/50;
		int cast_darea,count_darea=0;
		darea = (double*)calloc(10000, sizeof(double));
		p_darea = (double*)calloc(500000, sizeof(double));
		for (j=1; j<=10000; j++) {
			if(count3[j]==0 || j==CellNout[t-1])continue;
			//areal growth from previouos to current
			for (i=0; i<700; i++) {
				if(j==Cnumareacur1[i]){
					areag0=areacur1[i]/areaprev1[i];
				}
			}
			//areal growth from current to next
			for (i=0; i<700; i++) {
				if(j==Cnumareanext1[i]){
					areag1=areanext1[i]/areacur1[i];
				}
			}
			darea[j]=fabs(fabs(area1)-fabs(area0))/fabs(fabs(area1)+fabs(area0));
			//calculate probability density of darea
			if(isnan(darea[j])){
			}
			else{
				if(isinf(darea[j])){
				}
				else{
					bin1_darea=darea[j]/bin_darea;
					cast_darea=(int)bin1_darea;
					p_darea[cast_darea]=p_darea[cast_darea]+1.0;
					count_darea=count_darea+1;
				}
			}
		}
//Calculate 95% confidence interval
		double Bin_darea,sum_darea=0;
		for (j=0; j<50000; j++) {
			Bin_darea=(j+0.5)*bin_darea;
			sum_darea=sum_darea+p_darea[j]/(count_darea);
			if(sum_darea>0.95)thre_darea=Bin_darea;
			if(sum_darea>0.95)break;
		}
		for (j=1; j<10000; j++) {
			if(count3[j]==0 || j==CellNout[t-1])continue;
			for (k=0; k<700; k++) {
				if(j==parent[k]){ 
					for (i=0; i<N[t]; i++) {
						if(sun[k]==Cnum[i]){
							if(darea[j]<thre_darea)fprintf(file12,"%lf %lf %lf %lf\n",x1[i],y1[i],z1[i],darea[j]);
						}
					}
				}
			}
		}
		free(darea);
		free(p_darea);
		for (i=1; i<=16; i++)free(MTnear[i]);
		free(MTnear);
		free(x);
		free(y);
		free(z);
		free(Cnum);
		free(sun);
		free(parent);
		free(Sun);
		free(Parent);
		free(Cnumareaprev);
		free(areaprev);
		free(Cnumareacur);
		free(areacur);
		free(Cnumareanext);
		free(areanext);
		free(Sun3);
		free(Parent3);
		free(Cnumareaprev1);
		free(areaprev1);
		free(Cnumareacur1);
		free(areacur1);
		free(Cnumareanext1);
		free(areanext1);
	}
	
	return(0);
}
