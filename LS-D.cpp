
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>
#include <string.h>
#include <time.h>
#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>
#include <sstream>

#include <cstdlib>
#include "conio.h"
#include <time.h>
#include <fstream>
#include <iomanip>

using namespace std;

	int i,j,k,p;
	clock_t t1;


	int Nv,Np;
	int *IDp,*IDv;

	double bestcost=0,cost,cost1;
	double bestt;
	int *usedPM,*PM;
	int * VM;
	int*disk;
	int **N;
	double *C,*C2;
	int *sorted_c;
	int *tabp,*tabv;
	int voisinage=0;
	int decomp;int recherche;
	int count_neig;
	
	int * listeVm, *listePm;
	
	int *X,*X1,*Xbest;
	int **Y,**Y1,**Ybest;
	int *n_vm_par_pm,*n_vm_par_pm1,*n_vm_par_pm_best;
	double *res_CPU_phy,*res_CPU_phy1,*res_CPU_phy_best;
	double *res_MEM_phy,*res_MEM_phy1,*res_MEM_phy_best;
	double **res_DISK_phy,**res_DISK_phy1,**res_DISK_phy_best;
	
	
	int c[2];
	int v1,v2,p1,p2;
	string instance;
	int index;
	int cpt=0;
	int*y;
	float initcost;

	double CPU_phy[15]={8,8,8,8,16,16,16,16,16,32,48,64,80,120,120};
	double MEM_phy[15]={16,32,64,64,32,64,128,256,256,256,512,1024,2048,4096,4096};
	int nDISK_phy[15]={1,1,2,4,2,4,4,8,16,4,8,4,16,4,24};
	double COST_phy[15]={100,120,200,300,600,700,900,1500,1800,2500,3500,5000,7000,9000,1200};

	double DISK_phy[15][24]={{256,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
							{512,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
							{512,512,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
							{512,512,512,512,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
							{512,512,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
							{512,512,512,512,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
							{1000,1000,1000,1000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
							{1000,1000,1000,1000,1000,1000,1000,1000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
							{512,512,512,512,512,512,512,512,512,512,512,512,512,512,512,512,0,0,0,0,0,0,0,0},
							{1000,1000,1000,1000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
							{1000,1000,1000,1000,1000,1000,1000,1000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
							{1000,1000,1000,1000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
							{1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,0,0,0,0,0,0,0,0},
							{1000,1000,1000,1000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
							{1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600},
	};
							

	double CPU_vir[18]={1,2,4,8,2,4,8,16,32,2,4,8,16,32,4,8,16,32};
	double MEM_vir[18]={3.75,7.5,15,30,3.75,7.5,15,30,60,15.25,30.5,61,122,244,30.5,61,122,244};
	double REVENE_vir[18]={0.067,0.133,0.266,0.532,0.105,0.21,0.42,0.84,1.68,0.166,0.333,0.665,1.33,2.66,0.853,1.705,3.41,6.82};
	int nDISK_vir[18]={1,1,2,2,2,2,2,2,2,1,1,1,1,2,1,2,4,8};
	double DISK_vir[18][8]={
		{4,0,0,0,0,0,0,0},{32,0,0,0,0,0,0,0},{40,40,0,0,0,0,0,0},
		{80,80,0,0,0,0,0,0},{16,16,0,0,0,0,0,0},{40,40,0,0,0,0,0,0},
		{80,80,0,0,0,0,0,0},{160,160,0,0,0,0,0,0},{320,320,0,0,0,0,0,0},
		{32,0,0,0,0,0,0,0},{80,0,0,0,0,0,0,0},{160,0,0,0,0,0,0,0},
		{320,0,0,0,0,0,0,0},{320,320,0,0,0,0,0,0},{800,0,0,0,0,0,0,0},
		{800,800,0,0,0,0,0,0},{800,800,800,800,0,0,0,0},{800,800,800,800,800,800,800,800}
	};
	
	



double CPUp(int pm){	return CPU_phy[IDp[pm]];}
double COSTp(int pm){return COST_phy[IDp[pm]];}
double COSTUNITp(int pm){return (COST_phy[IDp[pm]]/CPUp(pm));}
double MEMp(int pm){	return MEM_phy[IDp[pm]];}
int nDISKp(int pm){	return nDISK_phy[IDp[pm]];}
double DISKp(int pm, int nd){ return DISK_phy[IDp[pm]][nd];}

double CPUv(int vm){	return CPU_vir[IDv[vm]];}
double MEMv(int vm){	return MEM_vir[IDv[vm]];}
double REVENUEv(int vm){return REVENE_vir[IDv[vm]];}
int nDISKv(int vm){	return nDISK_vir[IDv[vm]];}
double DISKv(int vm,int nd){ return DISK_vir[IDv[vm]][nd];}


double utilization(int pm){
	double res=1;
	double x1,x2,x3=0;
	x1=(res_CPU_phy[pm])/CPUp(pm);
	x2=(res_MEM_phy[pm])/MEMp(pm);
	for(k=0;k<nDISKp(pm);k++)
		x3+=res_DISK_phy[pm][k];
	x3=x3/(nDISKp(pm)*DISKp(pm,0));

	return x1*x2*x3*1000;
}

double utilization_add_vm(int vm,int pm){
	
	double x1,x2,x3=0;
	x1=(res_CPU_phy[pm]-CPUv(vm))/CPUp(pm);
	x2=(res_MEM_phy[pm]-MEMv(vm))/MEMp(pm);
	for(k=0;k<nDISKp(pm);k++)
		x3+=res_DISK_phy[pm][k];
	x3-=nDISKv(vm)*DISKv(vm,0);
	x3=(x3)/(nDISKp(pm)*DISKp(pm,0));

	return x1*x2*x3*1000;
}



void merge_VM(int,int,int);
void merge_sort_VM(int low,int high)
{
 int mid;
 if(low<high)
 {
  mid = low + (high-low)/2; //This avoids overflow when low, high are too large
  merge_sort_VM(low,mid);
  merge_sort_VM(mid+1,high);
  merge_VM(low,mid,high);
 }
}
void merge_VM(int low,int mid,int high)
{
 int h,i,j,*b,k;
 b=(int*) calloc(Nv,sizeof(int)); 

 h=low;
 i=low;
 j=mid+1;

 while((h<=mid)&&(j<=high))
 {
	 if(REVENUEv(VM[h])>=REVENUEv(VM[j]))
  {
   b[i]=VM[h];
   h++;
  }
  else
  {
   b[i]=VM[j];
   j++;
 }
  i++;
 }
 if(h>mid)
 {
  for(k=j;k<=high;k++)
  {
   b[i]=VM[k];
   i++;
  }
 }
 else
 {
  for(k=h;k<=mid;k++)
  {
   b[i]=VM[k];
   i++;
  }
 }
 for(k=low;k<=high;k++) VM[k]=b[k];
 free(b);
}

void merge_cost_PM(int low,int mid,int high)
{
 int h,i,j,*b,k;
 b=(int*) calloc(Np,sizeof(int)); 
 h=low;
 i=low;
 j=mid+1;

 while((h<=mid)&&(j<=high)) {
	 if(COSTp(PM[h])<=COSTp(PM[j])){
	   b[i]=PM[h];
	   h++;
	  }
	  else{
	   b[i]=PM[j];
	   j++;
	 }
  i++;
 }
 if(h>mid)
  for(k=j;k<=high;k++){
   b[i]=PM[k];
   i++;
  }
 else
  for(k=h;k<=mid;k++){
   b[i]=PM[k];
   i++;
  }
 
 for(k=low;k<=high;k++) 
	 PM[k]=b[k];
 free(b);
}
void merge_sort_cost_PM(int low,int high)
{
 int mid;
 if(low<high){
  mid = low + (high-low)/2; //This avoids overflow when low, high are too large
  merge_sort_cost_PM(low,mid);
  merge_sort_cost_PM(mid+1,high);
  merge_cost_PM(low,mid,high);
 }
}


void merge_PM(int low,int mid,int high)
{
 int h,i,j,k;
 int *b;
 b=(int*) calloc(Np,sizeof(int)); 
 h=low;
 i=low;
 j=mid+1;

 while((h<=mid)&&(j<=high)) {

if(COSTp(listePm[h])>=COSTp(listePm[j])){

	 
	   b[i]=listePm[h];
	   h++;
	  }
	  else{
	   b[i]=listePm[j];
	   j++;
	 }
  i++;
 }
 if(h>mid)
  for(k=j;k<=high;k++){
   b[i]=listePm[k];
   i++;
  }
 else
  for(k=h;k<=mid;k++){
   b[i]=listePm[k];
   i++;
  }
 
 for(k=low;k<=high;k++) 
	 listePm[k]=b[k];
 free(b);
}
void merge_sort_PM(int low,int high)
{
 int mid;
 if(low<high){
  mid = low + (high-low)/2; //This avoids overflow when low, high are too large
  merge_sort_PM(low,mid);
  merge_sort_PM(mid+1,high);
  merge_PM(low,mid,high);
 }
}


void merge_utilization_PM(int low,int mid,int high)
{
 int h,i,j,k;
 int *b;
 b=(int*) calloc(Np,sizeof(int)); 
 h=low;
 i=low;
 j=mid+1;

 while((h<=mid)&&(j<=high)) {
if(utilization(listePm[h])>=utilization(listePm[j])){
	   b[i]=listePm[h];
	   h++;
	  }
	  else{
	   b[i]=listePm[j];
	   j++;
	 }
  i++;
 }
 if(h>mid)
  for(k=j;k<=high;k++){
   b[i]=listePm[k];
   i++;
  }
 else
  for(k=h;k<=mid;k++){
   b[i]=listePm[k];
   i++;
  }
 
 for(k=low;k<=high;k++) 
	 listePm[k]=b[k];
 free(b);
}
void merge_sort_utilization_PM(int low,int high)
{
 int mid;
 if(low<high){
  mid = low + (high-low)/2; //This avoids overflow when low, high are too large
  merge_sort_utilization_PM(low,mid);
  merge_sort_utilization_PM(mid+1,high);
  merge_utilization_PM(low,mid,high);
 }
}


void merge_disk(int,int,int,int);
void merge_sort_disk(int low,int high,int pm)
{
 int mid;
 if(low<high)
 {
  mid = low + (high-low)/2; //This avoids overflow when low, high are too large
  merge_sort_disk(low,mid,pm);
  merge_sort_disk(mid+1,high,pm);
  merge_disk(low,mid,high,pm);
 }
}
void merge_disk(int low,int mid,int high,int pm)
{
 int h,i,j,*b,k;
 b=(int*) calloc(32,sizeof(int)); 

 h=low;
 i=low;
 j=mid+1;

 while((h<=mid)&&(j<=high))
 {
	 if(res_DISK_phy[pm][disk[h]]>=res_DISK_phy[pm][disk[j]])
  {
   b[i]=disk[h];
   h++;
  }
  else
  {
   b[i]=disk[j];
   j++;
 }
  i++;
 }
 if(h>mid)
 {
  for(k=j;k<=high;k++)
  {
   b[i]=disk[k];
   i++;
  }
 }
 else
 {
  for(k=h;k<=mid;k++)
  {
   b[i]=disk[k];
   i++;
  }
 }
 for(k=low;k<=high;k++) disk[k]=b[k];
 free(b);
}


int count_neighboor_structure(){

	int cn=decomp;
	int x=1;
	while(x<decomp){
		for(int cpt=0;cpt<decomp-x;cpt++)
			cn++;
		x++;
	}

	
	return cn;

}
void generate_table(){
	tabv=(int*) calloc ((decomp+1),sizeof(int));
	tabp=(int*) calloc ((decomp+1),sizeof(int));
	
	for(int cpt=0;cpt<decomp+1;cpt++){
		tabv[cpt]=(cpt*Nv)/decomp;
		tabp[cpt]=(cpt*Np)/decomp;

	}tabp[decomp]=Np-1;
	tabv[decomp]=Nv;

}
void init_recherche(){

	
	c[0]=1;c[1]=1;
	count_neig=count_neighboor_structure();
	generate_table();

	
	VM=(int*) calloc (Nv,sizeof(int));
	
	listeVm=(int*) calloc (Nv,sizeof(int));
	listePm=(int*) calloc (Np,sizeof(int));
	for(i=0;i<Np;i++) 
		listePm[i]=i;

	n_vm_par_pm=(int*) calloc (Np,sizeof(int));
	for(k=0;k<Nv;k++){
		n_vm_par_pm[X[k]]++;
	}

	Xbest=(int*) calloc (Nv,sizeof(int));

	Ybest=(int**) calloc (Nv,sizeof(int*));
	for(k=0;k<Nv;k++){
		Ybest[k]=(int*) calloc (nDISKv(k),sizeof(int));
	}
	

	for (int p=0;p<Nv;p++)
		Xbest[p]=X[p];

	for (int k=0;k<Nv;k++){
		for (int p=0;p<nDISKv(k);p++)
		Ybest[k][p]=Y[k][p];
	}

	n_vm_par_pm_best=(int*) calloc (Np,sizeof(int));
	for(int k=0;k<Np;k++){
		n_vm_par_pm_best[k]=n_vm_par_pm[k];
	}

	
	
	cost=bestcost;	
	res_CPU_phy_best=(double*) malloc (sizeof(double)*Np);
	res_MEM_phy_best=(double*) malloc (sizeof(double)*Np);
	res_DISK_phy_best=(double**) malloc (sizeof(double*)*Np);
	for(i=0;i<Np;i++)
		res_DISK_phy_best[i]=(double*) malloc (sizeof(double)*nDISKp(i));
	for(i=0;i<Np;i++){
		res_CPU_phy_best[i]=res_CPU_phy[i];
		res_MEM_phy_best[i]=res_MEM_phy[i];
		for(int k=0;k<nDISKp(i);k++)
			res_DISK_phy_best[i][k]=res_DISK_phy[i][k];
	}
	
}
void init_solution(){

	X=(int*) calloc (Nv,sizeof(int));

	Y=(int**) calloc (Nv,sizeof(int*));
	for(k=0;k<Nv;k++){
		Y[k]=(int*) calloc (8,sizeof(int));
	}
		
	res_CPU_phy=(double*) malloc (sizeof(double)*Np);
	res_MEM_phy=(double*) malloc (sizeof(double)*Np);
	res_DISK_phy=(double**) malloc (sizeof(double*)*Np);
	for(i=0;i<Np;i++)
		res_DISK_phy[i]=(double*) malloc (sizeof(double)*nDISKp(i));
	
	usedPM=(int*) calloc (Np,sizeof(int));
	PM=(int*) malloc (sizeof(int)*Np);
	VM=(int*) calloc (Nv,sizeof(int));
	for(i=0;i<Np;i++) 
		PM[i]=i;
	for(k=0;k<Nv;k++)
		VM[k]=k;
	


	
		
	for(i=0;i<Np;i++){
		res_CPU_phy[i]=CPUp(i);
		res_MEM_phy[i]=MEMp(i);
		for(k=0;k<nDISKp(i);k++)
			res_DISK_phy[i][k]=DISKp(i,k);
	}
	

}

void lecture_data_center(string fichier){
		
	fstream fdc(fichier, ios::in); 
	if(!fdc) printf("Error opening Data_center.txt file...\n");
		
	fdc >> Np;
	
	IDp=(int*) malloc (sizeof(int)*Np);

	
	for(i=0; i<Np;i++){
		fdc>>IDp[i];
		
	}

	fdc.close();
}
void lecture_requette(string fichier){
	
	fstream fr(fichier, ios::in);
	if(!fr) printf("Error opening %s.txt file...\n",fichier);
	
	fr>>Nv;

	IDv=(int*) malloc (sizeof(int)*Nv);

	for(i=0; i<Nv;i++){
		fr>>IDv[i];
	}

	fr.close();	
}

int* embedd_PM_VM(int i, int j){

	if((res_CPU_phy[i]<CPUv(j))||(res_MEM_phy[i]<MEMv(j))||(nDISKp(i)<nDISKv(j))) return nullptr;
	else{
		int *used_phys_disk=(int*) calloc (nDISKp(i),sizeof(int));
		disk=(int*) calloc (nDISKp(i),sizeof(int));
		for(k=0;k<nDISKp(i);k++) disk[k]=k;
		int *y=(int*) calloc (nDISKv(j),sizeof(int));
		

		int n=0;
		merge_sort_disk(0,nDISKp(i)-1,i);

		for(p=0;p<nDISKv(j);p++){
			
			for(k=0;k<nDISKp(i);k++){
				if((used_phys_disk[disk[k]]==0)&&(res_DISK_phy[i][disk[k]]>=DISKv(j,p))){
					y[p]=disk[k];
					used_phys_disk[disk[k]]=1;
					n++;
					
					break;
				}
			}
			
		}
		
		if(n<nDISKv(j)) { return nullptr;}
		else{
			res_CPU_phy[i]-=CPUv(j);
			res_MEM_phy[i]-=MEMv(j);
			X[j]=i;
			
			for(int p=0;p<nDISKv(j);p++){
				Y[j][p]=y[p];
				res_DISK_phy[i][y[p]]-=DISKv(j,p);
			}		
			free(disk);
			free(used_phys_disk);
			return y;
		}
	
	}
	
	return nullptr;
}

bool prim_test_embed_PM_VM(int i , int j){

	if((res_CPU_phy[i]<CPUv(j))||(res_MEM_phy[i]<MEMv(j))||(nDISKp(i)<nDISKv(j))){ return false;}
	else return true;
}

int* test_embedd_PM_VM(int i, int j){

	if((res_CPU_phy[i]<CPUv(j))||(res_MEM_phy[i]<MEMv(j))||(nDISKp(i)<nDISKv(j))){ return nullptr;}
	else{
		disk=(int*) calloc (nDISKp(i),sizeof(int));
		for(k=0;k<nDISKp(i);k++) disk[k]=k;
		int *used_phys_disk=(int*) calloc (nDISKp(i),sizeof(int));
		int *y=(int*) calloc (nDISKv(j),sizeof(int));
		

		int n=0;
		merge_sort_disk(0,nDISKp(i)-1,i);
		for(p=0;p<nDISKv(j);p++){
			for(k=0;k<nDISKp(i);k++){
				
				if((used_phys_disk[disk[k]]==0)&&(res_DISK_phy[i][disk[k]]>=DISKv(j,p))){
					y[p]=disk[k];
					used_phys_disk[disk[k]]=1;
					n++;
					
					break;
				}
			}
			
		}
		free(disk);
		free(used_phys_disk);

		if(n<nDISKv(j)) return nullptr;
		else return y;
		
	
	}
	
	return nullptr;
}

void remove_PM_VM(int i, int j, int*y){
	res_CPU_phy[i]+=CPUv(j);
	res_MEM_phy[i]+=MEMv(j);
	for(int p=0;p<nDISKv(j);p++)
		res_DISK_phy[i][y[p]]+=DISKv(j,p);	
		
}
void add_PM_VM(int i, int j, int*y){
		res_CPU_phy[i]-=CPUv(j);
		res_MEM_phy[i]-=MEMv(j);
		X[j]=i;
			
		for(int p=0;p<nDISKv(j);p++){
				Y[j][p]=y[p];
				res_DISK_phy[i][y[p]]-=DISKv(j,p);		
		}	
	
}

double calcul_cost(double c,int pm1,int pm2){
	//remove vm from pm1
	//add vm to pm2
	if(n_vm_par_pm[pm2]==0) c=c+COSTp(pm2);
	if(n_vm_par_pm[pm1]-1==0) c=c-COSTp(pm1);
	
	return c;
}
double calcul_cost_remove(double c,int pm){
	//remove vm from pm
	if(n_vm_par_pm[pm]-1==0) c=c-COSTp(pm);
	return c;
}
double calcul_cost_add(double c,int pm){
	//add vm to pm
	if(n_vm_par_pm[pm]==0) c=c+COSTp(pm);
	return c;
}
	
void calcul_bornes(int val,int r){
	int x=0;
		int va=val;
		while(va-decomp>=1){
			va=va-decomp;
			x++;
		}
	if(r==1){
		p1=tabp[decomp-x-1];p2=tabp[decomp-x];
		v1=tabv[va-1];v2=tabv[va];
	}else if(r==2){
		v1=tabv[decomp-va];v2=tabv[decomp-va+1];
		p1=tabp[decomp-va];p2=tabp[decomp-va+1];
	}else if(r==3){
		p1=tabp[decomp-va];p2=tabp[decomp-va+1];
	}
}


int nbR;
void recherche_local_N1(){

	cpt=0;
	
	for(i=0;i<Np;i++) 
		listePm[i]=i;
		
	merge_sort_PM(0,Np-1);
	
	for(i=0;i<Np;i++)
		for(k=0;k<Nv;k++)
			if(X[k]==listePm[i]){listeVm[cpt]=k;cpt++;}
			
	int index=0;
	double c1;
	int pm1;
	int vm1;
	
	for(int i=v1;i<v2;i++){
		index=0;
		for(int j=p2;j>=p1;j--){
			y=test_embedd_PM_VM(listePm[j],listeVm[i]);
			if((y!=nullptr)&&(listePm[j]!=X[listeVm[i]])){
				
				if((index==0) ||(utilization_add_vm(listeVm[i],listePm[j])-utilization(X[listeVm[i]])<c1)){
					c1=(utilization_add_vm(listeVm[i],listePm[j])-utilization(X[listeVm[i]]));
					pm1=j;
					vm1=i;
					index++;
				}

			}	
		}
			
		
		if(index!=0){

			cost=calcul_cost_add(cost,listePm[pm1]);
			cost=calcul_cost_remove(cost,X[listeVm[vm1]]);

			n_vm_par_pm[X[listeVm[vm1]]]--;
			n_vm_par_pm[listePm[pm1]]++;
		
			remove_PM_VM(X[listeVm[vm1]],listeVm[vm1],Y[listeVm[vm1]]);
			embedd_PM_VM(listePm[pm1],listeVm[vm1]);
		}
		
	}
}

void no_update_best(){
	cost=bestcost;
	for (int p=0;p<Nv;p++)
		X[p]=Xbest[p];

	for (int k=0;k<Nv;k++){
		for (int p=0;p<nDISKv(k);p++)
			Y[k][p]=Ybest[k][p];
	}
	for(int k=0;k<Np;k++)
		n_vm_par_pm[k]=n_vm_par_pm_best[k];
			
	for(int p=0;p<Np;p++){
		res_CPU_phy[p]=res_CPU_phy_best[p];
		res_MEM_phy[p]=res_MEM_phy_best[p];
		for(int k=0;k<nDISKp(p);k++)
			res_DISK_phy[p][k]=res_DISK_phy_best[p][k];
	}
			
	
}
void update_best(){
	
			
	printf("\n--------------UPDATE------------\n");	
			
	bestcost=cost;
	for (int p=0;p<Nv;p++)
		Xbest[p]=X[p];

	for (int k=0;k<Nv;k++)
		for (int p=0;p<nDISKv(k);p++)
			Ybest[k][p]=Y[k][p];
	for(int k=0;k<Np;k++)
		n_vm_par_pm_best[k]=n_vm_par_pm[k];
			
	for(int p=0;p<Np;p++){
		res_CPU_phy_best[p]=res_CPU_phy[p];
		res_MEM_phy_best[p]=res_MEM_phy[p];
		for(int k=0;k<nDISKp(p);k++)
			res_DISK_phy_best[p][k]=res_DISK_phy[p][k];
	}
			
}


bool FFD(){
	merge_sort_VM(0,Nv-1);

	bestcost=0;
	bool embed;

	for(j=0;j<Nv;j++){
		embed=false;
			
		//try to embedd vm on a used pm
		for(i=0;i<Np;i++){
			if (usedPM[i]==1)
				if(embedd_PM_VM(PM[i],VM[j])!=nullptr){
					embed=true;
					break;
				}
		}
		//try to embed vm on an unsed pm
		if(embed==false){
				
			for(i=0;i<Np;i++){
				if (usedPM[i]==0)
					if(embedd_PM_VM(PM[i],VM[j])!=nullptr){
						embed=true;
						bestcost=bestcost+COSTp(PM[i]);
						usedPM[i]=1;
						break;
					}
			}
		}
		if(embed==false){ printf("\n%d)NO embedd VM(%d) diskv:%d!",j, VM[j],nDISKv(VM[j]));break;}
		
	}
	return embed;
}


void main()
{

	fstream fl("instances_liste_large.txt", ios::in);
	if(!fl)	printf("Error opening file instances_liste.txt...\n");
	
	ofstream ffres("res.xls", ios::out | ios::app );
	ffres<<endl<<endl;
	
	
	while(!fl.eof()){
	fl>>instance;
	
	decomp=20;
	lecture_data_center("all_instances_large/data_center_"+instance);
	lecture_requette("all_instances_large/"+instance);
	

	t1 = clock();
	init_solution();
	printf("begin resolution...\n");

	FFD();
		
	double t = (double)(clock() - t1)/(double)CLK_TCK ;
	double t_best=t;
		
	printf("CPU=%f \ncost=%f\n\n\n",t,bestcost);
	ffres<<instance<<"\t"<<decomp<<"\t"<<bestcost<<"\t";
	
	
	init_recherche();		
	int iter=0;
	int val=1;

	double initcost=bestcost;
	

	recherche=1; val=1;
	do{
		calcul_bornes(val,recherche);
		
		recherche_local_N1();
		if(cost<bestcost) update_best();
		else no_update_best();
		val++;
	}while(val<=decomp*decomp);
		
	t = (double)(clock() - t1)/(double)CLK_TCK ;;	
	ffres<<bestcost<<"\t"<<t<<"\t"<<(initcost-bestcost)*100/initcost<<"\t";

	
	t = (double)(clock() - t1)/(double)CLK_TCK ;		
	printf("CPU=%0.3f \n\n",t);
	ffres<<bestcost<<"\t"<<t;
	ffres<<"\t"<<(initcost-bestcost)*100/initcost<<endl;


	merge_sort_cost_PM(0,Np-1);	
	
	for(i=0;i<Np;i++){
		if(n_vm_par_pm[i]>0) printf("%d ",i+1);
	}
	
	
	free(usedPM);free(VM);
	free(n_vm_par_pm);
	free(n_vm_par_pm_best);
	free(PM);
	free(tabp);free(tabv);
	free(X);free(Xbest);
	free(listeVm);free(listePm);
	free(res_CPU_phy);free(res_CPU_phy_best);
	free(res_MEM_phy);free(res_MEM_phy_best);

	for(i=0;i<Nv;i++){
		free(Y[i]);free(Ybest[i]);
	}
	free(Y);free(Ybest);
	for(i=0;i<Np;i++){
		free(res_DISK_phy[i]);free(res_DISK_phy_best[i]);
	}
	free(res_DISK_phy);free(res_DISK_phy_best);
	
	}
	

	ffres.close();

	system("pause");
	
	
}
