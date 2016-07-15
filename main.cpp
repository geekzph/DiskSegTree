//
//  main.cpp
//  multi-branch segment tree
//
//  Created by zph on 16/5/29.
//  Copyright © 2016年 zph. All rights reserved.
//
// not using precompoled header files

#include <iostream>
#include <ostream>
#include <fstream>
#include <stdio.h>
#include <vector>
#include <time.h>
#include <stdlib.h>
#include <malloc.h>
using namespace std;

const int kBranchNum = 10;                       //const branch number
fstream file;

typedef struct Info{
	int maxi, lmaxi, rmaxi, sum;
};
typedef struct Node{                             //Node structure
    int l;
    int r;
    int branch;
    int maxi,lmaxi,rmaxi,sum;
    int p[kBranchNum + 1];
    struct Node *s[kBranchNum + 1];
	Info info[kBranchNum][kBranchNum];
	int sequence;
} Node; 


typedef struct DiskNode{                             //Node structure
	int l;
	int r;
	int branch;
	int maxi, lmaxi, rmaxi, sum;
	int p[kBranchNum + 1];
} DiskNode;


typedef struct Kp{                              //pointer structer
    struct Node *s[kBranchNum+1];
} Kp;

typedef struct KpDisk{                              //pointer structer
	struct DiskNode *s[kBranchNum + 1];
} KpDisk;

Node *rootnode = (Node *)malloc(sizeof(Node)); //apply for root node
DiskNode *rootdisknode = (DiskNode *)malloc(sizeof(DiskNode)); //apply for root node
int g_node_size = sizeof(DiskNode);                //teh size of the Node structuer
int g_data_num = 0;                            //the amount of the dataset

float * ReadDate(char filename[])
{
    g_data_num = 0;
    float r[20004];
    FILE *fp;
    double result = 0.0;
    char strline[8];                            //each line's max-character
	errno_t err;
	err = fopen_s(&fp, filename, "r");
	if (err != 0)                               //open file
	{
		printf("error!");
		return NULL;
	}
    while (!feof(fp))
    {
        fgets(strline,12,fp);                   //read a line
        result=atof(strline);
        g_data_num++;
        r[g_data_num] =result;
        
        
    }
    fclose(fp);                                 //close file
    printf("read data file successed\n");
    return r;
}

//return max value
int max(int a, int b)
{
    return a>b?a:b;
}

//merge branch in CreateTree
void MergeBranchInC(Node *p,int n)
{
    int psum=0,pmaxi=0,plmaxi=0,prmaxi=0;
    for (int i = 1; i < n; i++)
    {
        if (i == 1)
        {
            p->sum = p->s[1]->sum + p->s[2]->sum;
            p->maxi = max(p->s[1]->maxi,max(p->s[2]->maxi,p->s[1]->rmaxi+p->s[2]->lmaxi));
            p->lmaxi = max(p->s[1]->lmaxi,p->s[1]->sum + p->s[2]->lmaxi);
            p->rmaxi = max(p->s[2]->rmaxi,p->s[2]->sum + p->s[1]->rmaxi);
            psum=p->sum ,pmaxi=p->maxi,plmaxi=p->lmaxi,prmaxi=p->rmaxi;
        }
        else if(i > 1)
        {
            int m = i + 1;
            p->sum = psum + p->s[m]->sum;
            p->maxi = max(pmaxi,max(p->s[m]->maxi,prmaxi+p->s[m]->lmaxi));
            p->lmaxi = max(plmaxi,psum + p->s[m]->lmaxi);
            p->rmaxi = max(p->s[m]->rmaxi,p->s[m]->sum + prmaxi);
            psum=p->sum ,pmaxi=p->maxi,plmaxi=p->lmaxi,prmaxi=p->rmaxi;
        }
        
    }
    
}

//merge branch in AddInfo
void MergeBranchInAddInfo(Node *p, int i, int j)
{
	int psum = 0, pmaxi = 0, plmaxi = 0, prmaxi = 0;
	int sum = 0,  maxi = 0, lmaxi = 0, rmaxi = 0, m = i;
	psum = p->s[i]->sum, pmaxi = p->s[i]->maxi, plmaxi = p->s[i]->lmaxi, prmaxi = p->s[i]->rmaxi;
	for (i = i ; i < j; i++)
	{
		int m = i + 1;
		sum = psum + p->s[m]->sum;
		maxi = max(pmaxi, max(p->s[m]->maxi, prmaxi + p->s[m]->lmaxi));
		lmaxi = max(plmaxi, psum + p->s[m]->lmaxi);
		rmaxi = max(p->s[m]->rmaxi, p->s[m]->sum + prmaxi);
		psum = sum, pmaxi = maxi, plmaxi = lmaxi, prmaxi = rmaxi;
	}
	p->info[m][j].sum = psum;
	p->info[m][j].maxi = pmaxi;
	p->info[m][j].lmaxi = plmaxi;
	p->info[m][j].rmaxi = prmaxi;
}

//merge branch in QuerySeg funcution
//QuerySeg funtion's merge procedure is different from CreateTree funtion merge procedure
Node *MergeBranchInQ(Node *p, Info *k, Node *q, int n)
{
    Node *res = (Node *)malloc(sizeof(Node));
        //only two branch to merge
    if(n ==1 )
    {
        res->sum = p->sum + q->sum;
        res->lmaxi = max(p->lmaxi,p->sum + q->lmaxi);
        res->rmaxi = max(q->rmaxi,q->sum + p->rmaxi);
        res->maxi = max(p->rmaxi+q->lmaxi,max(p->maxi,q->maxi));
        return res;

    }
	if (n == 3)
	{
		int psum = 0, pmaxi = 0, plmaxi = 0, prmaxi = 0;
		res->sum = p->sum + k->sum;
		res->lmaxi = max(p->lmaxi, p->sum + k->lmaxi);
		res->rmaxi = max(k->rmaxi, k->sum + p->rmaxi);
		res->maxi = max(p->rmaxi + k->lmaxi, max(p->maxi, k->maxi));

		psum = res->sum, pmaxi = res->maxi, plmaxi = res->lmaxi, prmaxi = res->rmaxi;

		res->sum = psum + q->sum;
		res->lmaxi = max(plmaxi, psum + q->lmaxi);
		res->rmaxi = max(q->rmaxi, q->sum + prmaxi);
		res->maxi = max(prmaxi + q->lmaxi, max(pmaxi, q->maxi));
	}
   
	return res;
}

//create multi-branch tree
void CreateTree(int l ,int r , Node *tp, float x[])
{
    int i = 1;
    int ll = l;
    int rr = r;
    tp->l = l;
    tp->r = r;
    int seg = (r - l) / kBranchNum;                   //segment length
    int flag = 0;
    
    if(ll == rr)
    {
        tp->l = ll;
        tp->r = rr;
        tp->maxi = tp->lmaxi = tp->rmaxi = tp->sum = x[ll];
		tp->branch = 0;
        return;
    }
    
    for (i = 1 ; i < kBranchNum + 1; i++)
    {
        seg = (r - l) / kBranchNum;
        if (ll > r || flag == 1)
        {
            break;
        }
        else if (ll + seg > r && ll < r + 1)
        {
            tp->s[i] = (Node *)malloc(sizeof(Node));   //apply for node
            CreateTree(ll, r,tp->s[i],x);              //create node
            flag = 1;                                  //quit for
        }
        else
        {
            tp->s[i] = (Node *)malloc(sizeof(Node));   //apply for node
            CreateTree(ll, ll + seg,tp->s[i],x);       //create node
            ll = ll + seg + 1;
            rr = ll + seg ;
        }
        tp -> branch = i;
    }
    MergeBranchInC(tp, i - 1);
}

void AddInfo(Node* root)
{
	if (!root)
	{
		return;
	}
	vector<Node*> vec;
	vec.push_back(root);

	int cur = 0;
	int last = 1;
	Node *p = root;
	while (cur<vec.size())
	{
		last = int(vec.size());
		while (cur<last)
		{
			p = vec[cur];
			cout << vec[cur]->l << "  ";
			//insert info 
			for (int i = 2; i <= vec[cur] -> branch - 1; i++)
			{
				for (int j = i; j <= vec[cur] -> branch - 1; j++)
				{
					MergeBranchInAddInfo(vec[cur], i, j);
					
				}
			}
			
			for (int i = 1; i <= p->branch; i++)
			{
				vec.push_back(vec[cur]->s[i]);
			}
			cur++;

		}
		cout << endl;
	}
}
//calculate pointer number
int IntervalNum(int x,Node *p)
{
    int m = 0;
    for(int i = 1;i <= p->branch; i++)
    {
        if (x >= p->s[i]->l && x <= p->s[i]->r) m = i;
    }
    return m;
}

//query segment from aa to bb's maxsub sum
Node *QuerySeg(int l,int r,int aa,int bb, Node *tp,int num)
{
    int flag1 = 0;                                    //only one branch
    int flag2 = 0;                                    //if should merge branch
    if(num != 0) tp = tp -> s[num];                   //change pointer
    Node *ka, *kl, *kr, *res, *res1;
    Kp *k = (Kp *)malloc(sizeof(Kp));                 //multi-branch's pointer
    res = (Node *)malloc(sizeof(Node));
    ka = (Node *)malloc(sizeof(Node));
    kl = (Node *)malloc(sizeof(Node));
    kr = (Node *)malloc(sizeof(Node));
	Info *info = new Info;
    if(aa <= l && bb >= r)
        return tp;
    int ll = 0;
    int rr = 0;
    ll = IntervalNum(aa, tp);                        //calculate pointer number
    rr = IntervalNum(bb, tp);
    if(ll == rr)
    {
        if(tp -> s[rr] -> r < bb)
            ka = QuerySeg(tp -> s[ll] ->r,tp -> s[rr]->r, aa, tp -> s[rr]->r, tp, ll);
        else
            ka = QuerySeg(tp -> s[ll] ->l,tp -> s[rr]->r, aa, bb, tp, ll);
        flag1 = 1;
    }
    int n = 1;
    if(ll < rr)
    {
        kl = QuerySeg(tp -> s[ll] ->l,tp -> s[ll] ->r, aa, tp -> s[ll] ->r, tp, ll);//lefmost point
        /*while (ll < rr - 1)
        {
            k -> s[i] = QuerySeg(tp -> s[ll+1] ->l,tp -> s[ll+1] ->r, tp -> s[ll+1] ->l, tp -> s[ll+1] ->r, tp, ll+1);
            ll++;
            i++;
        }*/
		if (ll < rr - 1)
		{
			info->maxi = tp->info[ll + 1][rr - 1].maxi;
			info->lmaxi = tp->info[ll + 1][rr - 1].lmaxi;
			info->rmaxi = tp->info[ll + 1][rr - 1].rmaxi;
			info->sum = tp->info[ll + 1][rr - 1].sum;
			n = 3;
		}
		
        kr = QuerySeg(tp -> s[rr] ->l,tp -> s[rr] ->r, tp -> s[rr] ->l, bb, tp, rr);//rightmost point
        flag2 = 1;
		
    }
    if (flag1 == 1)
    {
        res = ka;
    }
    else if(flag2 == 1)
    {
        res1 = MergeBranchInQ(kl, info, kr, n);
        res -> sum = res1 -> sum;
        res -> lmaxi = res1 -> lmaxi;
        res -> rmaxi = res1 -> rmaxi;
        res -> maxi = res1 -> maxi;
    }
    return res;
}

//merge branch in QuerySeg funcution
//QuerySeg funtion's merge procedure is different from CreateTree funtion merge procedure
DiskNode *MergeBranchInQDisk(DiskNode *p, KpDisk *k, DiskNode *q, int n)
{
	DiskNode *res = (DiskNode *)malloc(sizeof(DiskNode));
	//only two branch to merge
	if (n == 1)
	{
		res->sum = p->sum + q->sum;
		res->lmaxi = max(p->lmaxi, p->sum + q->lmaxi);
		res->rmaxi = max(q->rmaxi, q->sum + p->rmaxi);
		res->maxi = max(p->rmaxi + q->lmaxi, max(p->maxi, q->maxi));
		return res;

	}
	else
	{
		int psum = 0, pmaxi = 0, plmaxi = 0, prmaxi = 0;
		for (int i = 1; i < n; i++)
		{
			if (i == 1)
			{
				res->sum = p->sum + k->s[1]->sum;
				res->lmaxi = max(p->lmaxi, p->sum + k->s[1]->lmaxi);
				res->rmaxi = max(k->s[1]->rmaxi, k->s[1]->sum + p->rmaxi);
				res->maxi = max(p->rmaxi + k->s[1]->lmaxi, max(p->maxi, k->s[1]->maxi));
				psum = res->sum, pmaxi = res->maxi, plmaxi = res->lmaxi, prmaxi = res->rmaxi;
			}
			else if (i > 1)
			{
				res->sum = psum + k->s[i]->sum;
				res->lmaxi = max(plmaxi, psum + k->s[i]->lmaxi);
				res->rmaxi = max(k->s[i]->rmaxi, k->s[i]->sum + prmaxi);
				res->maxi = max(prmaxi + k->s[i]->lmaxi, max(pmaxi, k->s[i]->maxi));
				psum = res->sum, pmaxi = res->maxi, plmaxi = res->lmaxi, prmaxi = res->rmaxi;
			}

		}
		res->sum = psum + q->sum;
		res->lmaxi = max(plmaxi, psum + q->lmaxi);
		res->rmaxi = max(q->rmaxi, q->sum + prmaxi);
		res->maxi = max(prmaxi + q->lmaxi, max(pmaxi, q->maxi));
		return res;
	}

}
//calculate pointer number
DiskNode *p = new DiskNode;
int IntervalNumInDisk(int x, DiskNode *tp)
{
	
	int m = 0;
	for (int i = 1; i <= tp->branch; i++)
	{
		file.seekg(tp->p[i], ios::beg);     //change pointer
		file.read((char*)p, g_node_size);
		if (x >= p->l && x <= p->r) m = i;
	}
	//free(p);
	return m;
}



//query segment from aa to bb's maxsub sum in disk
DiskNode *QuerySegInDisk(int l, int r, int aa, int bb, DiskNode *tp, int num)
{
	DiskNode *tpl = new DiskNode;                      //tpl is the leftmost pointer
	DiskNode *tpr = new DiskNode;                      //tpm is the pointer between tpl and tpr
	DiskNode *tpm = new DiskNode;                      //tpr is the rightmost pointer
	DiskNode *tpa = new DiskNode;
	int flag1 = 0;                             //only one branch
	int flag2 = 0;                             //if should merge branch
	if (num != 0)
	{
		file.seekg(tp->p[num], ios::beg);      //change pointer
		file.read((char*)tpa, g_node_size);
	}
	
	else file.read((char*)tpa, g_node_size);
	DiskNode *ka = new DiskNode;
	DiskNode *kl = new DiskNode;
	DiskNode *kr = new DiskNode;
	DiskNode *res = new DiskNode;
	DiskNode *res1 = new DiskNode;
	KpDisk *k = new KpDisk;           //multi-branch's pointer
	if (aa <= l && bb >= r)
	{
		delete[] tpl;
		delete[] tpr;
		delete[] tpm;

		delete[] ka;
		delete[] kl;
		delete[] kr;
		delete[] res;
		delete[] res1;
		delete[] k;
		return tpa;
	}
		
	int ll = 0;
	int rr = 0;
	ll = IntervalNumInDisk(aa, tpa);            //calculate pointer number
	rr = IntervalNumInDisk(bb, tpa);
	file.seekg(tpa->p[ll], ios::beg);
	file.read((char*)tpl, g_node_size);
	file.seekg(tpa->p[rr], ios::beg);
	file.read((char*)tpr, g_node_size);
	if (ll == rr)
	{
		if (tpr->r < bb)
			ka = QuerySegInDisk(tpl->r, tpr->r, aa, tpr->r, tpa, ll);
		else
			ka = QuerySegInDisk(tpl->l, tpr->r, aa, bb, tpa, ll);
		flag1 = 1;
	}
	int i = 1;
	if (ll < rr)
	{
		kl = QuerySegInDisk(tpl->l, tpl->r, aa, tpl->r, tpa, ll);   //lefmost point
		while (ll < rr - 1)
		{
			file.seekg(tpa->p[ll + 1], ios::beg);
			file.read((char*)tpm, g_node_size);
			k->s[i] = QuerySegInDisk(tpm->l, tpm->r, tpm->l, tpm->r, tpa, ll + 1);
			ll++;
			i++;
		}
		kr = QuerySegInDisk(tpr->l, tpr->r, tpr->l, bb, tpa, rr);   //rightmost point
		flag2 = 1;

	}
	if (flag1 == 1)
	{
		res = ka;
	}
	else if (flag2 == 1)
	{
		res1 = MergeBranchInQDisk(kl, k, kr, i);
		res->sum = res1->sum;
		res->lmaxi = res1->lmaxi;
		res->rmaxi = res1->rmaxi;
		res->maxi = res1->maxi;
	}
	delete[] tpl;
	delete[] tpr;
	delete[] tpm;


	delete[] kl;
	delete[] kr;

	delete[] res1;
	return res;
}

//create index file
void WriteIndexFile(Node* root)
{
	fstream file;
	file.open("test.dat", ios::out | ios::trunc | ios::binary);  //create index file
	if (!file)
		cout << "error! can't create file" << endl;

    if (!root)
    {
        return;
    }
    vector<Node*> vec;
    vec.push_back(root);
    
    int cur=0;
    int last=1;
    Node *p = root;
    while (cur<vec.size())
    {
        last=int(vec.size());
        while (cur<last)
        {
			p = vec[cur];
            cout<<vec[cur]->l<<"  ";
            vec[cur] -> sequence = cur;
            for(int i = 1;i <= p->branch; i++)
            {
                vec.push_back(vec[cur]->s[i]);
            }
            cur++;
           
        }
        cout<<endl;
    }
    
    cur = 0;
    last = 1;
    p =root;
	DiskNode *disknode = new DiskNode;
    while (cur<vec.size())
    {
        last=int(vec.size());
        while (cur<last)
        {
			p = vec[cur];
			
            for(int i = 1;i <= p->branch; i++)
            {
                vec[cur] -> p[i] = vec[cur] -> s[i] -> sequence * sizeof(DiskNode);
            }
			disknode->l = vec[cur]->l;
			disknode->r = vec[cur]->r;
			disknode->branch = vec[cur]->branch;
			disknode->maxi = vec[cur]->maxi;
			disknode->lmaxi = vec[cur]->lmaxi;
			disknode->rmaxi= vec[cur]->rmaxi;
			disknode->sum = vec[cur]->sum;
			memcpy(disknode->p, vec[cur]->p, sizeof(int)*(kBranchNum + 1));
			//file.write((char*)vec[cur], sizeof(Node));
			file.write((char*)disknode, sizeof(DiskNode));
			file.seekp(0, ios::end);
            cur++;
            
        }
        cout<<endl;
    }
	file.seekp(0, ios::end);
	file.close();
}

//the function to test read from index file
void ReadNode()
{
	fstream file;
	file.open("test.dat", ios::in | ios::binary);           //open index file
	DiskNode *s = new DiskNode;
	file.read((char*)s, sizeof(DiskNode));
	cout << "zuobaian : " << s -> l << endl;
	cout << "youbian : " << s -> r << endl;
	cout << "pointer : " << s -> p[1] << endl;
	file.seekg(s -> p[1], ios::beg);
	file.read((char*)s, sizeof(DiskNode));
	cout << "zuobaian : " << s->l << endl;
	cout << "youbian : " << s->r << endl;
	cout << "pointer : " << s->p[1] << endl;
	file.close();
}

int main(int argc, const char * argv[]) {
    float *p;
    p = ReadDate("100.txt");                                //loda dataset
    printf("total data is %d\n",g_data_num);                //show the amount of the dataset
	int a = 1;                                              //a to b index
	int b = 100;
	int left = 1;
	int right = 100;
    CreateTree(a, b, rootnode, p);                        //create index in memory
	AddInfo(rootnode);
	//WriteIndexFile(rootnode);                             //write index to disk
	//ReadNode();
	file.open("test.dat", ios::in | ios::binary);           //open index file
	while (1)
	{
		cin >> left >> right;
		Node *res1 = QuerySeg(a, b, left, right, rootnode, 0);
		//DiskNode *res2 = QuerySegInDisk(a, b, left, right, rootdisknode, 0);//query from aa to bb 's max segment sum
		cout << "maxsub sum is " << res1->maxi << " in memory" << endl;
		//cout << "maxsub sum is " << res2->maxi << " in index file" << endl;
	}
	
	
	file.close();
    free(rootnode);
	free(rootdisknode);
    return 0;
}
