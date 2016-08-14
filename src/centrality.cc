/*
 * centrality.cc
 *
 *  Created on: 17 juin 2016
 *      Author: e-spin
 */

#include "centrality.h"

// Vector and find() library
#include "vector"
#include <algorithm>

// Cout library
#include <iostream>

// String stream dependencies:
#include <fstream>
#include <sstream>

//Set precision for cout
#include <iomanip>

//Random generator
#include <stdlib.h>

// Parallelization
#include <omp.h>
	int cpt=0;


void Network::BC(){
	int v;
	int i;
	int j;
	int w,thread_id;
	std::vector<int>::const_iterator it;
	double o=0;
	//	#pragma omp parallel for reduction(+:this->bc)
    // parallelize this chunk of code

	for(v=0;v<this->size;v++){
		std::vector<double> d(this->size,-1.0);
		std::vector<double> sigma(this->size,0.0);
		std::vector<int> S;
		std::vector<int> Q;
		std::vector<	std::vector<int> > P (this->size);
		sigma[v]=1;

		d[v]=0;
		Q.insert(Q.end(),v);
		while(!Q.empty()){
			j=Q[Q.size()-1];
			Q.pop_back();
			S.insert(S.end(),j);
			for(int id = 0; id !=this->getAdjsize(j); ++id){
				if(d[this->a[j][id]]<0){
					Q.insert(Q.begin(),this->a[j][id]);
					d[this->a[j][id]]=d[j]+1;
				}
				if(d[this->a[j][id]]==(d[j]+1)){
					sigma[this->a[j][id]]=sigma[this->a[j][id]]+sigma[j];
					P[this->a[j][id]].insert(P[this->a[j][id]].end(),j);
				}
			}
		}
		std::vector<double> delta(this->size,0.0);
		while(!S.empty()){
			w=S[S.size()-1];
			S.pop_back();

			for(it= P[w].begin();it!= P[w].end(); ++it) {
				delta[*it]=delta[*it]+sigma[*it]/sigma[w]*(1+delta[w]);

			}
			if(w!=v){
				this->bc[w]=this->bc[w]+delta[w]/this->size;
			}
		}
		if(double(v)/double(this->size)>=o){
			std::cout<<"\n Process en cours : "<<100*o<<" % done !";
			o=o+0.1;
		}
	}

//	for(i=0;i<this->size;i++){
//		this->bc[i]=this->bc[i]/this->size;
//	}

}

void Network::BC_E(){
	int v;
	int i;
	int j;
	int w,thread_id;

	std::vector<int>::const_iterator it;
	double o=0;
	//	#pragma omp parallel for reduction(+:this->bc)
    // parallelize this chunk of code

	for(v=0;v<this->size;v++){
		std::vector<double> d(this->size,-1.0);
		std::vector<double> sigma(this->size,0.0);
		std::vector<int> S;
		std::vector<int> Q;
		std::vector<	std::vector<int> > P (this->size);
		sigma[v]=1;

		d[v]=0;
		Q.insert(Q.end(),v);
		while(!Q.empty()){
			j=Q[Q.size()-1];
			Q.pop_back();
			S.insert(S.end(),j);
			for(int id = 0; id !=this->getAdjsize(j); ++id){
				if(d[this->a[j][id]]<0){
					Q.insert(Q.begin(),this->a[j][id]);
					d[this->a[j][id]]=d[j]+1;
				}
				if(d[this->a[j][id]]==(d[j]+1)){
					sigma[this->a[j][id]]=sigma[this->a[j][id]]+sigma[j];
					P[this->a[j][id]].insert(P[this->a[j][id]].end(),j);
				}
			}
		}
		std::vector<double> delta1(this->size,0.0);
		std::vector<double> delta2(this->size,0.0);
		std::vector<double> bc1(this->size,0);
		std::vector<double> bc2(this->size,0);

		while(!S.empty()){
			w=S[S.size()-1];
			S.pop_back();

			for(it= P[w].begin();it!= P[w].end(); ++it) {
				delta1[*it]=delta1[*it]+sigma[*it]/sigma[w]*(d[w]*sigma[w]+delta1[w]);
				delta2[*it]=delta2[*it]+sigma[*it]/sigma[w]*(sigma[w]+delta2[w]);
			}
//			if(w!=v){
				bc1[w]=bc1[w]+delta1[w];
				bc2[w]=bc2[w]+delta2[w];
	//			if(delta2[w]!=0)
	//				this->bc[w]=this->bc[w]+delta1[w]/delta2[w];
//			}
		}
		for(i=0;i<this->size;i++){
			if(bc2[i]!=0)
				this->bc[i]=bc[i]+bc1[i]/bc2[i]/this->size;

		}
		if(double(v)/double(this->size)>=o){
			std::cout<<"\n Process en cours : "<<100*o<<" % done !";
			o=o+0.1;
		}
	}
}


void Network::BC_pond(){
	int v;
	std::pair<int,double> j;
	int w;
	std::vector<int>::const_iterator it;
	double o=0;

	for(v=0;v<this->size;v++){
		int id;
		std::vector<double> d(this->size,-1);
		std::vector<double> sigma(this->size,0.0);
		std::vector<int> S;
		std::vector<std::pair<int,double> > Q;
		std::vector<std::vector<int> > P(this->size);
		sigma[v]=1;
		double we;
//		NodeWeight[v]=0;

		d[v]=0;
		Q.insert(Q.begin(),std::pair<int,double>(v,0));
		int min;
		while(!Q.empty()){
			min = getMinPair(Q);
			j=Q[min];
			Q.erase(Q.begin()+min);
			S.insert(S.end(),j.first);

			for(id = 0; id !=this->getAdjsize(j.first); ++id){
				we=d[j.first]+weight[j.first][id]; //dist(x[a[j.first][id]],y[a[j.first][id]],x[j.first],y[j.first]);
				if(f(we-d[a[j.first][id]])<pow(10,-14)){
					sigma[a[j.first][id]]+=sigma[j.first];
					P[this->a[j.first][id]].insert(P[this->a[j.first][id]].end(),j.first);
				}
				if(sigma[a[j.first][id]]==0){
					Q.insert(Q.begin(),std::pair<int,double>(a[j.first][id],we));
				}
				else
					if(we<d[a[j.first][id]]){
						sigma[a[j.first][id]]=0;
						Q[getFirstPair(Q,a[j.first][id])].second=we;
						P[a[j.first][id]].clear();
					}else
						if(we>d[a[j.first][id]])
							continue;

				d[a[j.first][id]]=we;
				sigma[a[j.first][id]]+=sigma[j.first];
				P[this->a[j.first][id]].insert(P[this->a[j.first][id]].end(),j.first);
			}
		}
		std::vector<double> delta1(this->size,0.0);
		std::vector<double> delta2(this->size,0.0);
		std::vector<double> bc1(this->size,0);
		std::vector<double> bc2(this->size,0);

		while(!S.empty()){
			w=S[S.size()-1];
			S.pop_back();

			for(it= P[w].begin();it!= P[w].end(); ++it) {
				delta1[*it]=delta1[*it]+sigma[*it]/sigma[w]*(d[w]*sigma[w]+delta1[w]);
				delta2[*it]=delta2[*it]+sigma[*it]/sigma[w]*(sigma[w]+delta2[w]);
			}
//			if(w!=v){
				bc1[w]=bc1[w]+delta1[w];
				bc2[w]=bc2[w]+delta2[w];
	//			if(delta2[w]!=0)
	//				this->bc[w]=this->bc[w]+delta1[w]/delta2[w];
//			}
		}
		for(int i=0;i<this->size;i++){
			if(bc2[i]!=0)
				this->bc[i]=bc[i]+bc1[i]/bc2[i]/this->size;

		}
		/*
		std::vector<double> delta(this->size,0.0);
		while(!S.empty()){
			w=S[S.size()-1];
			S.pop_back();

			for(it= P[w].begin();it!= P[w].end(); ++it) {
				delta[*it]=delta[*it]+sigma[*it]/sigma[w]*(1+delta[w]);
			}
			if(w!=v){
				this->bc[w]=this->bc[w]+delta[w]/this->size;
			}
		}
		*/
//		if(double(v)/double(this->size)>=o){
//			std::cout<<"\n Process en cours : "<<100*o<<" % done !";
//			o=o+0.1;
//		}
	}
}

void Network::Percolation(double p){
	int i,j,size2, o;

	for(i=0;i<this->size;i++){
		size2=this->getAdjsize(i);
		o=0;
		for(j=0;j<size2;j++){
			if(((double) rand() / (RAND_MAX))<p){
				this->erase(this->a[i][o],i);
				o=o-1;
			}
			o=o+1;
		}
	}

}

void Network::GiantComponent(){
	std::vector<int> cx(this->size,0);
	int i,v,j,k,cpt,max,giant_component;

	giant_component=0;
	max=0;

	std::vector<int> Q;
	k=0;

	for(i=0;i<cx.size();i++){
		v=i;
		if(cx[v]==0){
			cpt=1;
			k=k+1;
			Q.insert(Q.end(),v);
			cx[v]=k;
			while(!Q.empty()){
				v=Q[Q.size()-1];
				Q.pop_back();
				for(j=0;j<this->getAdjsize(v);j++){
					if(cx[this->a[v][j]]==0){
						cpt+=1;
						Q.insert(Q.end(),this->a[v][j]);
						cx[this->a[v][j]]=cx[v];
					}

				}
			}
			if(cpt>max){
				max=cpt;
				giant_component=k;
			}
		}
	}
	int size, o;
	for(i=0;i<this->size;i++){
		o=0;
		if(cx[i]!=giant_component){
			size=this->getAdjsize(i);
			for(j=0;j<size;j++){
				this->erase(this->a[i][o],i);
			}
		}
	}
}

void Network::FilterBC(double u){
	int i,j,s;
	int cpt=0;
	this->NbNode=0;

	for(i=0;i<this->size;i++){
		if(this->bc[i]<u){
			this->bc[i]=0;
			cpt+=1;
			s=this->getAdjsize(i);
			for(j=0;j<s;j++){
				this->erase(this->a[i][0],i);
			}
		}
		else{
			this->bc[i]=1;
			this->NbNode+=1;
		}

	}

}

void Network::FilterBC(double u1,double u2){
	int i,j,s;
	int cpt=0;
	this->NbNode=0;

	for(i=0;i<this->size;i++){
		if(this->bc[i]<u1){
			this->bc[i]=0;
			cpt+=1;
			s=this->getAdjsize(i);
			for(j=0;j<s;j++){
				this->erase(this->a[i][0],i);
			}
		}
		else{
			this->bc[i]=1;
			this->NbNode+=1;
		}

	}

}

void Network::SweepLine(){
	int i,j;

	std::vector<size_t> idx(this->getXsize());
	std::vector<double> X;
	std::vector<double> Y;
	std::vector<std::vector<int> > A2(this->size);
	std::vector<double> nm;
	const std::vector<double> &v =x;

	std::vector<std::pair<int,int> > line;

	A2=this->a;
	for (size_t i = 0; i != idx.size(); ++i){
		idx[i] = i;
	}

//	sort(idx.begin(), idx.end(),[&X2](size_t i1, size_t i2) {return (sqrt(pow(X2[i1][0]-X2[i2][0],2)+(pow(X2[i2][1]-X2[i1][1],2)))<pow(10,-8));});
	sort(idx.begin(), idx.end(),[&v](size_t i1, size_t i2) {return (v[i1]<v[i2]);});

	X=this->x;
	Y=this->y;

	for(i=0;i<this->getXsize();i++){

		this->x[i]=X[idx[i]];
		this->y[i]=Y[idx[i]];
		for(j=0;j<this->getAdjsize(i);j++){
			this->a[i][j]=find(idx.begin(),idx.end(),A2[i][j])-idx.begin();
		}
	}
	A2=this->a;
	for(i=0;i<idx.size();i++){
		this->a[i]=A2[idx[i]];
	}
//	line.insert(line.begin(),std::pair<int,int>(i,this->a[i]));


	int k;
	int p1; int p3; int p2; int p4;
	int si=this->size;
/*
	for(i=0;i<s;i++){
		for(j=0;j<this->getAdjsize(i);j++){
			if(a[i][j]>i){
				line.insert(line.end(),std::pair<int,int>(i,this->a[i][j]));
				for(k=0;k<line.size()-1;k++){
					if(CheckPlanarity2(line[k].first,line[k].second,i,this->a[i][j])==1){
						int tmp;
						tmp=line[k].first;
						line[k].first=this->size-1;
						line[line.size()-1].second=this->size-1;
					}
				}
			}

			else{
				std::vector<std::pair<int,int> >::iterator it;
				int &v=a[i][j];
				int &vv=i;
				it=std::find_if( line.begin(), line.end(), [&v,&vv](const std::pair<int,int>& element){
					return (element.first==v)&&(element.second==vv);});
				for(k=s;k<this->size;k++){
					for(int w=0; w<this->getAdjsize(k);w++){
						CheckPlanarity2((*it).first,(*it).second,k,this->a[k][w]);
					}
				}
				// && a[i][j]<i
				if(it!=line.end()){
					line.erase(it);
				}
			}
*/
	std::vector<int> q;
	std::vector<std::vector <int> > intersection;

	for(i=0;i<si;i++){
		for(j=0;j<this->getAdjsize(i);j++){
			if(a[i][j]>i){
				line.insert(line.end(),std::pair<int,int>(i,this->a[i][j]));
				for(k=0;k<line.size()-1;k++){
					q=CheckPlanarity2(line[k].first,line[k].second,i,this->a[i][j]);
					if(q.size()!=0){
						intersection.insert(intersection.begin(),q);
					}
				}
			}

			else{
				std::vector<std::pair<int,int> >::iterator it;
				int &v=a[i][j];
				int &vv=i;
				it=std::find_if( line.begin(), line.end(), [&v,&vv](const std::pair<int,int>& element){
					return (element.first==v)&&(element.second==vv);});
				if(it!=line.end()){
					line.erase(it);
				}
			}
		}
	}


	for(i=0;i<intersection.size();i++){
		p1=intersection[i][3];
		p2=intersection[i][2];
		p3=intersection[i][1];
		p4=intersection[i][0];

		double xd,yd,xa,xb,xc,ya,yb,yc;

		xa=this->x[p1];
		ya=this->y[p1];

		xb=this->x[p2];
		yb=this->y[p2];

		xc=this->x[p3];
		yc=this->y[p3];

		xd=this->x[p4];
		yd=this->y[p4];



		double s1_x, s1_y, s2_x, s2_y;

		s1_x=xb-xa;     s1_y=yb-ya;
		s2_x=xd-xc;     s2_y=yd-yc;

		double s, t;
		s=(-s1_y *(xa-xc)+s1_x*(ya-yc))/(-s2_x*s1_y+s1_x*s2_y);
		t=(s2_x*(ya-yc)-s2_y*(xa-xc))/(-s2_x*s1_y+s1_x*s2_y);

		a.insert(a.end(),intersection[i]);
		this->size++;
		this->x.insert(x.end(),x[p1]+t*(x[p2]-x[p1]));
		this->y.insert(y.end(),y[p1]+t*(y[p2]-y[p1]));

		for(j=0;j<2;j++){
			p1=intersection[i][2*j+1];
			p2=intersection[i][2*j];
			if(find(a[p1].begin(),a[p1].end(),p2) != a[p1].end() && find(a[p2].begin(),a[p2].end(),p1) != a[p2].end()){
				a[p1][find(a[p1].begin(),a[p1].end(),p2)-a[p1].begin()]=a.size()-1;
				a[p2][find(a[p2].begin(),a[p2].end(),p1)-a[p2].begin()]=a.size()-1;
			}
			else{
				double dR=dist(x[p1],y[p1],x[this->size-1],y[this->size-1]);
				double dL=dist(x[p2],y[p2],x[this->size-1],y[this->size-1]);
				int nR=p1;
				int nL=p2;
				for(k=si;k<this->size-1;k++){
					if(AreAligned(p2,k,this->size-1)){
						if(dist(x[k],y[k],x[this->size-1],y[this->size-1])<dL){
							dL=dist(x[k],y[k],x[this->size-1],y[this->size-1]);
							nL=k;
						}
					}
					if(AreAligned(p1,k,this->size-1)){
						if(dist(x[k],y[k],x[this->size-1],y[this->size-1])<dR){
							dR=dist(x[k],y[k],x[this->size-1],y[this->size-1]);
							nR=k;
						}
					}
				}
				a[nR][find(a[nR].begin(),a[nR].end(),nL)-a[nR].begin()]=this->size-1;
				a[nL][find(a[nL].begin(),a[nL].end(),nR)-a[nL].begin()]=this->size-1;
				a[this->size-1][find(a[this->size-1].begin(),a[this->size-1].end(),p1)-a[this->size-1].begin()]=nR;
				a[this->size-1][find(a[this->size-1].begin(),a[this->size-1].end(),p2)-a[this->size-1].begin()]=nL;

			}
		}

	}
}


std::vector<int> Network::CheckPlanarity2(int p1, int p2, int p3, int p4){

	/*
	 * 	A(p1) linked vith B(p2)
	 *
	 * 	C(p3) linked with D(p4)
	 *
	 * 	Check if A-B and C-D cross :
	 *
	 */

	std::vector <int> q;

	int i,j,k,l;
	double xd,yd,xa,xb,xc,ya,yb,yc;

	xa=this->x[p1];
	ya=this->y[p1];

	xb=this->x[p2];
	yb=this->y[p2];

	xc=this->x[p3];
	yc=this->y[p3];

	xd=this->x[p4];
	yd=this->y[p4];



	double s1_x, s1_y, s2_x, s2_y;

	s1_x=xb-xa;     s1_y=yb-ya;
	s2_x=xd-xc;     s2_y=yd-yc;

	double s, t;
	s=(-s1_y *(xa-xc)+s1_x*(ya-yc))/(-s2_x*s1_y+s1_x*s2_y);
	t=(s2_x*(ya-yc)-s2_y*(xa-xc))/(-s2_x*s1_y+s1_x*s2_y);


	if ((s>pow(10,-5)) && (s<0.9999999) && (t>pow(10,-5)) && (t<0.9999999) && (p2!=p4)){
		cpt++;
		std::cout<<"Intersection"<<std::endl;
		std::cout<<std::fixed << std::setprecision(14) <<s<<"\t"<<t<<std::endl;

		q.clear();
		q.insert(q.begin(),p1);
		q.insert(q.begin(),p2);
		q.insert(q.begin(),p3);
		q.insert(q.begin(),p4);
		std::cout<<p1<<"\t"<<p2<<"\t"<<p3<<"\t"<<p4<<std::endl;
		return q;
	}
	return q;
}

void Network::CheckPlanarity(){

	/*
	 * 	A linked vith D
	 *
	 * 	B linked with C
	 *
	 */

	int i,j,k,l;

	std::vector<size_t> idx(this->getXsize());
	std::vector<double> X;
	std::vector<double> Y;
	std::vector<std::vector<int> > A2(this->size);
	std::vector<double> nm;
	const std::vector<double> &v =x;

	std::vector<int> mark(this->getNbEdge());

	std::vector<std::pair<int,int> > line;

	A2=this->a;
	for (size_t i = 0; i != idx.size(); ++i){
		idx[i] = i;
	}

//	sort(idx.begin(), idx.end(),[&X2](size_t i1, size_t i2) {return (sqrt(pow(X2[i1][0]-X2[i2][0],2)+(pow(X2[i2][1]-X2[i1][1],2)))<pow(10,-8));});
	sort(idx.begin(), idx.end(),[&v](size_t i1, size_t i2) {return (v[i1]<v[i2]);});

	X=this->x;
	Y=this->y;

	for(i=0;i<this->getXsize();i++){

		this->x[i]=X[idx[i]];
		this->y[i]=Y[idx[i]];
		for(j=0;j<this->getAdjsize(i);j++){
			this->a[i][j]=find(idx.begin(),idx.end(),A2[i][j])-idx.begin();
		}
	}
	A2=this->a;
	for(i=0;i<idx.size();i++){
		this->a[i]=A2[idx[i]];
	}

	getNodeList();
	double xd,yd,xa,xb,xc,ya,yb,yc;
	cpt=0;

	for(i=0;i<this->NodeList.size();i++){
		for(k=0;k<this->getAdjsize(this->NodeList[i]);k++){
			for(j=i+1;j<this->NodeList.size();j++){
				for(l=0;l<this->getAdjsize(this->NodeList[j]);l++){
					xa=this->x[NodeList[i]];
					ya=this->y[NodeList[i]];

					xb=this->x[NodeList[j]];
					yb=this->y[NodeList[j]];

					xd=this->x[this->a[this->NodeList[i]][k]];
					yd=this->y[this->a[this->NodeList[i]][k]];

					xc=this->x[this->a[this->NodeList[j]][l]];
					yc=this->y[this->a[this->NodeList[j]][l]];

				    double s1_x, s1_y, s2_x, s2_y;

				    s1_x=xd-xa;     s1_y=yd-ya;
				    s2_x=xc-xb;     s2_y=yc-yb;

				    double s, t;
				    s=(-s1_y *(xa-xb)+s1_x*(ya-yb))/(-s2_x*s1_y+s1_x*s2_y);
				    t=(s2_x*(ya-yb)-s2_y*(xa-xb))/(-s2_x*s1_y+s1_x*s2_y);


				    if ((s>pow(10,-5)) && (s<0.9999999) && (t>pow(10,-5)) && (t<0.9999999) && (this->a[this->NodeList[i]][k]!=this->a[this->NodeList[j]][l])){
						cpt++;
				    	std::cout<<"Intersection"<<std::endl;
						std::cout<<std::fixed << std::setprecision(14) <<s<<"\t"<<t<<std::endl;
						std::vector <int> q;
						int node_a,node_b,node_c,node_d;

						node_a=NodeList[i];
						node_b=NodeList[j];
						node_d=this->a[this->NodeList[i]][k];
						node_c=this->a[this->NodeList[j]][l];
						q.clear();
						q.insert(q.begin(),node_a);
						q.insert(q.begin(),node_b);
						q.insert(q.begin(),node_c);
						q.insert(q.begin(),node_d);
						std::cout<<node_a<<"\t"<<node_d<<"\t"<<node_b<<"\t"<<node_c<<std::endl;
						a[a.size()-1]=q;
						this->x.insert(x.end(),xa+t*(xd-xa));
						this->y.insert(y.end(),ya+t*(yd-ya));

						a[node_a][find(a[node_a].begin(),a[node_a].end(),node_d)-a[node_a].begin()]=a.size()-1;
						a[node_d][find(a[node_d].begin(),a[node_d].end(),node_a)-a[node_d].begin()]=a.size()-1;
						a[node_b][find(a[node_b].begin(),a[node_b].end(),node_c)-a[node_b].begin()]=a.size()-1;
						a[node_c][find(a[node_c].begin(),a[node_c].end(),node_b)-a[node_c].begin()]=a.size()-1;
						this->size++;
				    }
				}
			}

		}
	}
	std::cout<<"Compteur : "<<cpt<<std::endl;
}

void Network::DoubleEdge(){
	for(int i=0;i<a.size();i++){
		for(int j=0;j<a[i].size();j++){
			for(int k=0;k<count(a[i].begin(),a[i].end(),a[i][j])-1;k++){
				a[i].erase(find(a[i].begin(),a[i].end(),a[i][j]));
				std::cout<<"ERASE"<<std::endl;
			}
		}
	}
}



void Network::SortArray(){
	int i,j;

	std::vector<size_t> idx(this->getXsize());
	std::vector<double> X;
	std::vector<double> Y;
	std::vector<std::vector<int> > A2(this->size);
	std::vector<double> nm;
	const std::vector<double> &v =x;
	A2=this->a;

	for (size_t i = 0; i != idx.size(); ++i){
		idx[i] = i;
	}

	sort(idx.begin(), idx.end(),[&v](size_t i1, size_t i2) {return (v[i1]<v[i2]);});

	int s;
	X=this->x;
	Y=this->y;

	for(i=0;i<this->getXsize();i++){

		this->x[i]=X[idx[i]];
		this->y[i]=Y[idx[i]];
		for(j=0;j<this->getAdjsize(i);j++){
			this->a[i][j]=find(idx.begin(),idx.end(),A2[i][j])-idx.begin();
		}
	}
	A2=this->a;
	for(i=0;i<idx.size();i++){
		this->a[i]=A2[idx[i]];
	}

	int w;
	this->size=this->a.size();
	i=0;
	while(i<this->size){
		j=i+1;
		while(x[j]==x[i]){
			if(y[i]==y[j]){
				for(w=0;w<this->getAdjsize(j);w++){
					this->a[i].insert(this->a[i].begin(),this->a[j][w]);
					this->a[this->a[j][w]].erase(find(this->a[this->a[j][w]].begin(),this->a[this->a[j][w]].end(),j));
					this->a[this->a[j][w]].insert(this->a[this->a[j][w]].begin(),i);
				}
				a[j].clear();
			}
			j=j+1;
		}
		i=i+1;
	}
}

void Network::starNetwork(double Nb,double Nn){
	this->a.clear();
	this->x.clear();
	this->y.clear();
	this->NBranches=Nb;
	this->Nn=Nn;
	this->size=Nb*Nn+1;
	double X,Y;
	int i,j;

	create_node(Nn,Nn);

	for (i=0;i<Nn;i++){
		for(j=0;j<Nb;j++){

//			X=(Nn-1)*cos(j*M_PI/(Nb-1))+Nn*sin(j*M_PI/(Nb-1));
//			Y=-(Nn-1)*sin(j*M_PI/(Nb-1))+(Nn)*cos(j*M_PI/(Nb-1));

			X=-(i+1)*cos(j*2*M_PI/(Nb))+Nn;
			Y=-(i+1)*sin(j*2*M_PI/(Nb))+Nn;
//			Y=(int)Y+1*((Y-(int)Y)>0);

			create_node(X,Y);
			if(i>0){
				this->link(a.size()-1,a.size()-1-Nb,1);
//				a[a.size()-1].insert(a[a.size()-1].begin(),a.size()-1-Nb);
//				a[a.size()-1-Nb].insert(a[a.size()-1-Nb].begin(),a.size()-1);
			}
			else{
				this->link(a.size()-1,0,1);
			}
		}
		for(j=0;j<i;j++){

		}
	}
	this->set_bc();
}

void Network::setBoucle(int R, double w){
	/*
	 * R should be between 0 and Nn
	 */

	int i;
	double X,Y;
//	a[NBranches*(R+1)].insert(a[NBranches*(R+1)].begin(),NBranches*R+1);
//	a[NBranches*R+1].insert(a[NBranches*R+1].begin(),NBranches*(R+1));

	link(NBranches*(R-1)+1,NBranches*R,w);

	for(i=0;i<this->NBranches-1;i++){
//		a[i+NBranches*R+1].insert(a[i+NBranches*R+1].begin(),i+NBranches*R);
//		a[i+NBranches*R].insert(a[i+NBranches*R].begin(),i+NBranches*R+1);

		link(i+NBranches*(R-1)+1,i+NBranches*(R-1)+2,w);
	}

/*
	X=-R*cos(M_PI/(NBranches))+Nn;
	Y=R*sin(M_PI/(NBranches))+Nn;
	create_node(X,Y);
	link(NBranches*(R-1)+1,a.size()-1,w);
	link(a.size()-1,NBranches*R,w);
	this->size++;
	for(i=0;i<this->NBranches-1;i++){
		X=-R*cos((2*i+1)*M_PI/(NBranches))+Nn;
		Y=-R*sin((2*i+1)*M_PI/(NBranches))+Nn;
		create_node(X,Y);
		this->size++;
		link(i+NBranches*(R-1)+1,a.size()-1,w);
		link(a.size()-1,i+NBranches*(R-1)+2,w);
	}
*/

}



void Network::autoset(int K){
	this->size=K*K;
	double L=K;
	double X;
	int i;
	for (i=0;i<this->size;i++){
		X=i;
		create_node(fmod(X,L),(X-fmod(X,L))/L);
		if(i>=K){
			this->a[i].insert(this->a[i].end(),i-K);
		}
		if((i%K)>0){
			this->a[i].insert(this->a[i].end(),i-1);
		}
		if((i%K)<(K-1)){
			this->a[i].insert(this->a[i].end(),i+1);
		}
		if(i<(this->size-K)){
			this->a[i].insert(this->a[i].end(),i+K);
		}
	}
	this->set_bc();
}


void Network::set_adj(const char * string){
	this->a.clear();
    int i;
    std::stringstream temp_str;
    temp_str<<"../data/"<<this->name<<"/"<<string;
    std::string str = temp_str.str();
    const char* cstr2 = str.c_str();

    std::ifstream in(cstr2,std::ios::in);
    i=0;
    std::string line;
    std::vector<int> nm;
    while (std::getline(in,line)) {
            std::stringstream ls(line );
            int number;
            this->a.insert(this->a.end(),nm);
            while (ls>> number ) {
                    this->a[i].insert(this->a[i].begin(),number);
            }
            i=i+1;
    }
    std::cout<<i;
    this->size=i;
}

void Network::set_coord(const char * string){
	this->x.clear();
	this->y.clear();

    std::stringstream temp_str;
    temp_str<<"../data/"<<this->name<<"/"<<string;
    std::string str = temp_str.str();
    const char* cstr2 = str.c_str();

    std::ifstream in(cstr2,std::ios::in);

    std::string line;
    double x,y;
    while (in >> x>>y){
            this->x.insert(this->x.end(),x);
            this->y.insert(this->y.end(),y);
    }
}

void Network::set_bc(const char * string){
	this->bc.clear();

    std::stringstream temp_str;
    temp_str<<"../data/"<<this->name<<"/"<<string;
    std::string str = temp_str.str();
    const char* cstr2 = str.c_str();

    std::ifstream in(cstr2,std::ios::in);

    std::string line;
    double x,y,b;
    while (in >>x>>y>>b){
            this->bc.insert(this->bc.end(),b);
    }
}



void Network::plot_net(int w){

	int o;
	std::stringstream temp_str;
	temp_str<<"../data/"<<this->name<<"/net_"<<w<<".dat";
	std::string str = temp_str.str();
	const char* cstr2 = str.c_str();

	std::ofstream myfile(cstr2);
	if(myfile.is_open()){
		int j;

		std::vector<std::vector<int> >::iterator itr;
		o=0;

		for(itr=this->a.begin();itr!=this->a.end();itr++){
			for(j=0;j<this->getAdjsize(o);j++){
//				double x=X[id[o]];
//				double y=Y[id[o]];
				double x=this->x[o];
				double y=this->y[o];
				myfile<<std::fixed << std::setprecision(14) <<x<<"\t"<<y<<std::endl;
				x=this->x[(*itr)[j]];
				y=this->y[(*itr)[j]];
				myfile<<std::fixed << std::setprecision(14)<<x<<"\t"<<y<<std::endl;
				myfile<<std::endl;
			}
			o=o+1;
		}
	}
    myfile.close();

}

void Network::plot_periph(int O){
	std::stringstream temp_str;
	temp_str<<"../data/"<<this->name<<"/periph_"<<O<<".dat";
	std::string str = temp_str.str();
	const char* cstr2 = str.c_str();

	std::ofstream myfile(cstr2);
	if(myfile.is_open()){
		int i;
		for(i=1;i<this->periph.size();i++){
			myfile<<std::fixed << std::setprecision(14) <<this->x[this->periph[i-1]]<<"\t"<<this->y[this->periph[i-1]]<<std::endl;
			myfile<<std::fixed << std::setprecision(14) <<this->x[this->periph[i]]<<"\t"<<this->y[this->periph[i]]<<std::endl<<std::endl;

		}
	}
}
void Network::plot_coord(int w){

	std::stringstream temp_str;
	temp_str<<"../data/"<<this->name<<"/coord_"<<w<<".dat";
	std::string str = temp_str.str();
	const char* cstr2 = str.c_str();

	std::ofstream myfile(cstr2);
	if(myfile.is_open()){

		int i;
		for(i=0;i<this->getXsize();i++){
			double x=this->x[i];
			double y=this->y[i];
			myfile<<std::fixed << std::setprecision(14) <<x<<"\t"<<y<<std::endl;
		}
	}
    myfile.close();

}


void Network::plot_adj(int w){

	int o;
	std::stringstream temp_str;
	temp_str<<"../data/"<<this->name<<"/adj_"<<w<<".dat";
	std::string str = temp_str.str();
	const char* cstr2 = str.c_str();

	std::ofstream myfile(cstr2);
	if(myfile.is_open()){
		int j;

		std::vector<std::vector<int> >::iterator itr;
		o=0;

		for(itr=this->a.begin();itr!=this->a.end();itr++){
			for(j=0;j<this->getAdjsize(o);j++){
				myfile<<(*itr)[j]<<"\t";
			}
			myfile<<std::endl;
			o=o+1;
		}
	}
    myfile.close();

}

void Network::plot_bc(int w){

	int o;
	std::stringstream temp_str;
	temp_str<<"../data/"<<this->name<<"/bc_"<<w<<".dat";
	std::string str = temp_str.str();
	const char* cstr2 = str.c_str();

	std::ofstream myfile(cstr2);
	if(myfile.is_open()){
		std::vector<std::vector<int> >::iterator itr;
		o=0;

		if(w==-1){
			for(itr=this->a.begin();itr!=this->a.end();itr++){
				double x=this->x[o];
				double y=this->y[o];
				double b=this->bc[o];
				if((*itr).size()!=0){
					myfile<<std::fixed << std::setprecision(14) <<x<<"\t"<<y<<"\t"<<b<<std::endl;
					myfile<<std::endl;
				}
				else{
					myfile<<std::fixed << std::setprecision(14) <<x<<"\t"<<y<<"\t"<<0<<std::endl;
					myfile<<std::endl;
				}
				o=o+1;
			}

		}
		else{
			for(o=0;o<bc.size();o++){
				double x=this->x[o];
				double y=this->y[o];
				double b=this->bc[o];
				if(b!=0){
					myfile<<std::fixed << std::setprecision(14) <<x<<"\t"<<y<<"\t"<<b<<std::endl;
					myfile<<std::endl;
				}
			}

		}
	}
    myfile.close();

}

void Network::plot_vect(std::vector<double> vect, const char * string){
	std::stringstream temp_str;
	temp_str<<"../data/"<<this->name<<"/"<<string<<".dat";
	std::string str = temp_str.str();
	const char* cstr2 = str.c_str();

	std::ofstream myfile(cstr2);
	if(myfile.is_open()){
		int i;
		for(i=0;i<vect.size();i++){
			myfile<<i*this->step<<"\t"<<vect[i]<<"\n";
		}
	}
}

void Network::plot_bcPond(std::vector<double> vect, const char * string,double pond){
	std::stringstream temp_str;
	temp_str<<"../data/"<<this->name<<"/"<<string<<"_"<<pond<<".dat";
	std::string str = temp_str.str();
	const char* cstr2 = str.c_str();

	std::ofstream myfile(cstr2);
	if(myfile.is_open()){
		int i;
		int j=1;
		for(i=0;i<Nn+1;i++){
//			myfile<<j<<"\t"<<i%Nn+1<<"\t"<<vect[i]<<"\n";
//			if(i%Nn==(Nn-1)){
//				j=j+1;
//				myfile<<"\n";
//			}
			myfile<<i<<"\t"<<vect[i*NBranches]<<"\n";

		}
	}
}


void Network::plot_measures(){
	int i;

	for(i=this->loop-1;i>=0;i--){
		this->val[i]=this->val[i]/this->val[0];
		if (L[i]!=0)
				this->form[i]=this->area[i]/(M_PI*(0.5*this->L[i])*(0.5*this->L[i]));
//		std::cout<<L[i]<<"\t"<<area[i];
		this->area[i]=this->area[i]/this->area[0];
	}

	/*
	for(i=this->loop-1;i>=0;i--){
		this->L[i]=M_PI*(0.5*this->L[i])*(0.5*this->L[i]);
	}
	*/
	this->plot_vect(this->L,"area_norm");
	this->plot_vect(this->val,"perimeter");
	this->plot_vect(this->area,"area");
	this->plot_vect(this->form,"form");
}

void Network::filter(double x_min,double x_max, double y_min, double y_max ){
	int i,j;
	int o=0;
	for(i=0;i<this->size;i++){
		std::cout<<o<<std::endl;
		if(this->x[i]<x_min || this->x[i]>x_max || this->y[i]<y_min || this->y[i]>y_max){
				for(j=0;j<this->getAdjsize(i);j++){
						this->a[this->a[i][j]].erase(find(this->a[this->a[i][j]].begin(),this->a[this->a[i][j]].end(),i));
				}
				this->a[i].clear();
		}
		o++;
	}
}

void Network::RemoveEmptyNodes(){
		int i,j;
		this->getNodeList();
		int o=0;

		std::vector<double> X(NodeList.size());
		std::vector<double> Y(NodeList.size());
		std::vector<std::vector<int> > A(NodeList.size());
		std::vector<double> Bc(NodeList.size());
		std::vector<double> nm;

		for(i=0;i<NodeList.size();i++){
			X[i]=x[NodeList[i]];
			Y[i]=y[NodeList[i]];
			A[i]=a[NodeList[i]];
			if(this->bc.size()!=0)
				Bc[i]=bc[NodeList[i]];
			for(j=0;j<this->getAdjsize(NodeList[i]);j++){
				A[i][j]=find(NodeList.begin(),NodeList.end(),a[NodeList[i]][j])-NodeList.begin();
			}
		}
		if(this->bc.size()!=0){
			bc.clear();
			this->bc=Bc;
		}
		a.clear();
		x.clear();
		y.clear();
		this->a=A;
		this->x=X;
		this->y=Y;

/*
		for(i=0;i<this->size;i++){
			if(this->getAdjsize(i)==0 && i<this->NodeList.size()){
				o=o+1;
				int NodeSwitch=this->NodeList[this->NodeList.size()-o];
				this->a[i]=a[NodeSwitch];

				for(j=0;j<this->getAdjsize(i);j++){
					this->a[this->a[i][j]][find(a[a[i][j]].begin(),this->a[this->a[i][j]].end(),NodeSwitch)-a[a[i][j]].begin()]=i;
				}

				this->a.erase(a.begin()+NodeSwitch);

				this->x[i]=this->x[NodeSwitch];
				this->x.erase(x.begin()+NodeSwitch);

				this->y[i]=this->y[NodeSwitch];
				this->y.erase(y.begin()+NodeSwitch);

				this->bc[i]=this->bc[NodeSwitch];
				this->bc.erase(bc.begin()+NodeSwitch);
			}
		}
		i=a.size()-1;
		while(a[i].size()==0){
			a.erase(a.begin()+i);
			this->x.erase(x.begin()+i);
			this->y.erase(y.begin()+i);
			this->bc.erase(bc.begin()+i);
			i--;
		}
		this->size=a.size();
		*/
		this->size=a.size();
}

/*
 * 	MEASURES :
 */

void Network::initMeasures(int l){
	this->val=std::vector<double>(l,0);
	this->form=std::vector<double>(l,0);
	this->area=std::vector<double>(l,0);
	this->L=std::vector<double>(l,0);
	this->degree=std::vector<int>(this->size,0);
}

void Network::RemBranches(){
	int i,a,b;
	int k=1;
	int v,w;
	int res;
	std::vector<int> mark;
	std::vector<int> rem;

	std::vector<std::vector<int> >::iterator ite;
	this->degree.clear();
	for(ite=this->a.begin();ite!=this->a.end();ite++){
		this->degree.insert(this->degree.end(),(*ite).size());
		if((*ite).size()==0)
			mark.insert(mark.end(),1);
		else
			mark.insert(mark.end(),0);
	}
	this->CleanNet();
	k=2;

	for(int v=0;v<mark.size();v++){
		if(mark[v]==0){
			mark[v]=1;
			int s=this->getAdjsize(v);
			int p=0;
			for(i=0;i<s;i++){
				if(mark[this->a[v][p]]==0){
					w=1;
					a=v;
					b=this->a[v][p];
					this->erase(v,b);
					res=this->ShortestPath(a,b);
					if(res!=0){
						this->a[v].insert(this->a[v].end(),b);
						this->a[b].insert(this->a[b].end(),v);
					}
					else
						this->CleanNet();

				}
			}
		}
	}
}

double Network::maxLength(){

	int j,i;

	this->maxX=this->maxY=0;
	this->minX=this->minY=this->size;
	double tmp,max;
	max=0;

	for(j=0;j<this->size;j++){
			if(this->getAdjsize(j)>0){
				if(this->x[j]>maxX)
					maxX=this->x[j];
				if(this->y[j]>maxY)
					maxY=this->y[j];
				if(this->x[j]<minX)
					minX=this->x[j];
				if(this->y[j]<minY)
					minY=this->y[j];
//				for(i=j;i<this->size;i++){
	//				tmp=this->dist(x[i],y[i],x[j],y[j]);
		//			if(tmp>max)
			//			max=tmp;
			//	}
			}
		}
	max=sqrt((maxX-minX)*(maxX-minX)+(maxY-minY)*(maxY-minY));


	return max;

	/*
		if((maxX-minX)>(maxY-minY))
			return(maxX-minX);
		else
			return(maxY-minY);
	*/
}

int Network::ShortestPath(int a, int b){
	std::vector<int>length(this->size,-1);
	int i,o,j,cpt;
	cpt=0;
	std::vector<int> Q;
	std::vector<int> P;


	for(i=0;i<this->getAdjsize(a);i++){
		Q.insert(Q.end(),this->a[a][i]);
		length[this->a[a][i]]=1;
		if(this->a[a][i]==b){
			break;
		}
	}
	Nb=this->size;
	o=1;
	while(length[b]==-1 && cpt<Nb){
		cpt=cpt+1;
		o=o+1;
		for(i=0;i<Q.size();i++){
			for(j=0;j<this->getAdjsize(Q[i]);j++){
				if(length[this->a[Q[i]][j]]==-1){
					P.insert(P.end(),this->a[Q[i]][j]);
					length[this->a[Q[i]][j]]=o+1;
				}
			}
		}
		Q=P;
		P.clear();
	}
	if(length[b]!=-1){
		return 1;
	}
	return 0;
}

void Network::CleanNet(){
	int i,b,c,j,o;

	for(i=0;i<this->size;i++){
		if(this->degree[i]==1){
			j=this->a[i][0];
			o=i;
			this->erase(j,o);
			this->degree[o]=-1;

			this->Nb=this->Nb-1;

			this->degree[j]=this->degree[j]-1;

			while(this->degree[j]==1){
				o=j;
				j=this->a[o][0];
				this->erase(j,o);
				this->degree[o]=-1;
				this->Nb=this->Nb-1;
				this->degree[j]=this->degree[j]-1;
			}
		}
	}
}

void Network::Periph(int O){
	int i,mx,left,right,current;
	left=right=mx=0;
	double tmp,mn;
	mn=0;
	double max=0;
	double min=this->size;
	double amn=2*M_PI;
	double amx=0;

	for(i=0;i<this->size;i++){
		if(this->y[i]>max && this->getAdjsize(i)!=0){
			max=this->y[i];
			mx=i;
		}
		if(this->y[i]<min && this->getAdjsize(i)!=0){
			min=this->y[i];
			mn=i;
		}
	}

	std::cout<<mn<<"\t"<<mx<<std::endl;
	for(i=0;i<this->getAdjsize(mx);i++){
		tmp=this->angle(this->x[mx],this->x[this->a[mx][i]],this->y[mx],this->y[this->a[mx][i]]);
		if(tmp<amn){
			amn=tmp;
//	Left
			right=this->a[mx][i];
		}
		if(tmp>amx){
			amx=tmp;
			left=this->a[mx][i];
		}
	}
	current=mx;
	this->a[left].erase(find(this->a[left].begin(),this->a[left].end(),current));
	this->a[right].erase(find(this->a[right].begin(),this->a[right].end(),current));
	this->a[current].clear();
	this->periph.insert(this->periph.begin(),current);
	int o=0;
	mn=angle(this->x[left],this->x[current],this->y[left],this->y[current]);
	while(left!=right ){
		current=left;
		amn=2*M_PI;
		for(i=0;i<this->getAdjsize(current);i++){
			tmp=angle(this->x[current],this->x[this->a[current][i]],this->y[current],this->y[this->a[current][i]])-mn;
			if((tmp+2*M_PI*(tmp<0))<amn){
				amn=tmp+2*M_PI*(tmp<0);
				left=this->a[current][i];
			}
		}

		mn=angle(this->x[left],this->x[current],this->y[left],this->y[current]);
		this->a[left].erase(find(this->a[left].begin(),this->a[left].end(),current));
		this->a[current].erase((find(this->a[current].begin(),this->a[current].end(),left)));
		this->periph.insert(this->periph.begin(),current);

		o=o+1;
//		std::cout<<o<<std::endl;
//		this->plot_periph(0);
	}

	this->periph.insert(periph.begin(),right);
	this->periph.insert(periph.begin(),mx);

	this->val[O]=0;
	for(i=1;i<this->periph.size();i++){
		this->val[O]=this->val[O]+dist(this->x[this->periph[i-1]],this->y[this->periph[i-1]],this->x[this->periph[i]],this->y[this->periph[i]]);
//		this->area[O]=this->area[O]+(this->x[this->periph[i]]-this->x[this->periph[i-1]])*this->f(this->y[this->periph[i]]-this->y[this->periph[i-1]]);
		this->area[O]=this->area[O]+0.5*(this->x[this->periph[i]]-this->x[this->periph[i-1]])*(this->y[this->periph[i]]+this->y[this->periph[i-1]]);
	}
}

