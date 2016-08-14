/*
 * centrality.h
 *
 *  Created on: 17 juin 2016
 *      Author: e-spin
 */

#ifndef SRC_CENTRALITY_H_
#define SRC_CENTRALITY_H_

#include <../vector>
#include <../math.h>
#include <algorithm>

// String stream dependencies:
#include <fstream>
#include <sstream>

// STD
#include <iostream>

class Network
{

public:

/*
 * Network variables:
 */

	int NBranches;
	int Nn;
	double R;
	std::vector<double> bc;
	std::vector<int> degree;
	std::vector<double> x;
	std::vector<double> y;
	std::vector< std::vector <int> > a;
	std::vector< std::vector <double> > weight;
	double step;
	int loop;
	int size;
	int Nb;
	int NbEdge;
	int NbNode;
	std::vector<int> NodeList;
	const char * name;

	double maxX, maxY, minX, minY;

/*
 * Measures variables:
 */

	std::vector<double> periph;
	std::vector<double> val;
	std::vector<double> area;
	std::vector<double> L;
	std::vector<double> form;

/*
 * Initialisation:
 */
	void starNetwork(double,double);
	void setBoucle(int,double);
	void autoset(int);
	void initMeasures(int);

	void set_adj(const char * );
	void set_coord(const char * );
	void set_bc(const char * );
	void set_bc(){
		this->bc.clear();
		this->bc = std::vector<double> (this->size);
	}
	void set_name(const char * n){this->name=n;};

/*
 * Main network functions:
 */

	void DoubleEdge();
	void RemoveEmptyNodes();
	void ThreadBc();
	void BC_E();
	void BC_pond();
	void BC();
	void FilterBC(double);
	void FilterBC(double,double);
	void GiantComponent();
	void Percolation(double);
	void SortArray();
	void filter(double,double,double,double);

/*
 * Measures function:
 */
	void SweepLine();
	void CheckPlanarity();
	std::vector<int> CheckPlanarity2(int,int,int,int);
	void Periph(int);
	void RemBranches();
	void CleanNet();
	int ShortestPath(int,int);

/*
 * 	Plot functions:
 */

	void plot_adj(int);
	void plot_coord(int);
	void plot_bc(int);
	void plot_bcPond(std::vector<double>, const char *,double pond);
	void plot_net(int);
	void plot_periph(int );
	void plot_measures();
	void plot_vect(std::vector<double>, const char *);
/*
 * Small functions:
 */

	void create_node(double Xnode, double Ynode){
		this->x.insert(this->x.end(),Xnode);
		this->y.insert(this->y.end(),Ynode);
		this->a.insert(this->a.end(),std::vector <int> ());
		this->weight.insert(this->weight.end(),std::vector <double> ());
	}

	int getFirstPair(std::vector< std::pair<int,double> > vect,int item){
		int o;
		o=0;
		while(vect[o].first!=item && o!=(vect.size()-1)){
			o=o+1;
		}
		return o;
	}


	void erase(int index,int item){
		this->a[index].erase(find(this->a[index].begin(),this->a[index].end(),item));
		this->a[item].erase(find(this->a[item].begin(),this->a[item].end(),index));
	}

	void erase_w(int index,int item){
		int j=find(this->a[index].begin(),this->a[index].end(),item)-a[index].begin()-1;
		this->weight[index].erase(weight[index].begin()+j);
		j=find(this->a[item].begin(),this->a[item].end(),index)-a[item].begin()-1;
		this->weight[item].erase(weight[item].begin()+j);
		this->a[index].erase(find(this->a[index].begin(),this->a[index].end(),item));
		this->a[item].erase(find(this->a[item].begin(),this->a[item].end(),index));
	}

	int getItem(std::vector<int> vect,int item){
		int size=vect.size();
		int o;
		o=0;
		while(o!=size){
			if(vect[o]==item)
				return o;
			else o+=1;
		}
		return -1;
	}

	int getMinPair(std::vector< std::pair<int,double> > vect){
		int size=vect.size();
		int o;
		int index;
		double min;
		o=0;
		min=10000;
		while(o!=size){
			if(vect[o].second<min){
				min=vect[o].second;
				index=o;
			}
			o=o+1;
		}
		return index;
	}

	double maxElement(std::vector<double> b){
		int i;
		int s=b.size();
		double m=0;

		for(i=0;i<s;i++){
			if(b[i]>m)
				m=b[i];
		}
		return m;
	}

	double maxElementIndex(std::vector<double> b){
		int i,k;
		int s=b.size();
		double m=0;
		k=0;

		for(i=0;i<s;i++){
			if(b[i]>m){
				m=b[i];
				k=i;
			}
		}
		return k;
	}

	double maxLength();

	double dist(double x1,double y1,double x2, double y2){
		return(sqrt(pow(x1-x2,2)+pow(y1-y2,2)));
	}

	void link(int p1,int p2, double w){
		a[p1].insert(a[p1].begin(),p2);
		a[p2].insert(a[p2].begin(),p1);

		if(w==0){
			weight[p1].insert(weight[p1].begin(),dist(x[p1],y[p1],x[p2],y[p2]));
			weight[p2].insert(weight[p2].begin(),dist(x[p1],y[p1],x[p2],y[p2]));
		}
		else{
			weight[p1].insert(weight[p1].begin(),w);
			weight[p2].insert(weight[p2].begin(),w);
		}
	}

	double angle(double xa,double xb,double ya,double yb){
		double angle;
		double X,Y;
		X=xb-xa;
		Y=yb-ya;
		if(X>0){
			if(Y>=0)
				angle=atan(Y/X);
			else
				angle=3*M_PI/2+f(atan(X/Y));
		}
		else{
			if(Y>0)
				angle=M_PI/2+f(atan(X/Y));
			else
				angle=M_PI+f(atan(Y/X));
		}
		return angle;
	}

	double projection(double xa, double ya, double xb, double yb, double xc, double yc){
		double t,m1,m2;
		double length;
		double X,Y,xd,yd;

		if(xc==xb){
			length=ya-yb;
			return length;
		}

		if(yc==yb){
			length=xa-xc;
			return length;
		}

		if(xc!=xb && yc!=yb){
		t=(yc-yb)/(xc-xb);

		m1=ya+1/t*xa;
		m2=(yb-t*xb);

		xd=(1/(1/t+t))*(m1-m2);
		yd=t*xd+(yb-t*xb);

		length=(xa-xd-(ya-yd)*t)*(1/sqrt(1+t*t));
		return length;
		}
		return 0;
	}
	const char * tostring(const char * string, int i){
		std::stringstream temp_str;
		temp_str<<string<<i;
		std::string str = temp_str.str();
		return str.c_str();
	}

	double f(double x){
		return x*(1-2*(x<0));
	}

	int getXsize(){
		return this->x.size();
	}
	int getYsize(){
		return this->y.size();
	}
	int getAdjsize(){
		return this->a.size();
	}
	int getAdjsize(int i){
		return this->a[i].size();
	}
	int getCbsize(){
		return this->bc.size();
	}
	int getNbNode(){
		int i;
		this->NbNode=0;
		for(i=0;i<this->a.size();i++)
			if(this->a[i].size()!=0){
				NbNode+=1;
			}
		return this->NbNode;
	}
	int getNbEdge(){
		int i;
		this->NbEdge=0;
		for(i=0;i<this->a.size();i++)
			NbEdge+=this->a.size();
		return this->NbEdge;
	}


//	template<typename T, typename... Args>
//	T sum(T a, Args... args) { return a + sum(args...); }

	void getNodeList(){
		int i;
		this->NodeList.clear();
		for(i=0;i<this->a.size();i++){
			if(this->a[i].size()!=0){
				NodeList.insert(NodeList.end(),i);
			}
		}
	}

    int checkDirectionality(){
            int i,j;
            for(i=0;i<this->a.size();i++){
                    for(j=0;j<this->a[i].size();j++){
                            if(find(a[a[i][j]].begin(),a[a[i][j]].end(),i)==a[a[i][j]].end()){
                                    std::cout<<i<<std::endl;
                                    a[a[i][j]].insert(a[a[i][j]].end(),i);
                            }
                    }
            }
            return 0;
    }

    void plot(std::vector<double> vect, const char * string){
    	std::stringstream temp_str;
    	temp_str<<"../data/"<<this->name<<"/"<<string<<".dat";
    	std::string str = temp_str.str();
    	const char* cstr2 = str.c_str();

    	std::ofstream myfile(cstr2);
    	if(myfile.is_open()){
    		int i;
    		for(i=0;i<vect.size();i++){
    			myfile<<i<<"\t"<<vect[i]<<"\n";
    		}
    	}
    }
    void removeNodes(std::vector<double> vect){

    	for(int i=0;i<vect.size();i++){
    		for(int j=0;j<a[vect[i]].size();j++){
    			a[a[vect[i]][j]].erase(find(a[a[vect[i]][j]].begin(),a[a[vect[i]][j]].end(),vect[i]));
    		}
    		a[vect[i]].clear();
    	}

    }

    int AreAligned(int p1,int p2,int p3){

    	int i;
    	double s;
    //	int p[3];
    	std::vector<double> d;
    	double epsilon=pow(10,-14);

    	d.insert(d.end(),dist(x[p1],y[p1],x[p3],y[p3]));
    	d.insert(d.end(),dist(x[p2],y[p2],x[p3],y[p3]));
    	d.insert(d.end(),dist(x[p2],y[p2],x[p1],y[p1]));

  //  	i=maxElementIndex(d);
  //  	p1=p[i%3];
 //   	p2=p[(i+1)%3];
 //   	p3=p[(i+2)%3];

//    	s=d[0]+d[1]+d[2];
		if( f(d[0]-(d[1]+d[2]))<epsilon){
			return 1;
		}
		return 0;
    }

};




#endif /* SRC_CENTRALITY_H_ */
