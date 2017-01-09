//#include<tools.h>

using namespace std;



template<class T> void MSmergeLists(T* list1, T* list2, T* target, int len1, int len2, bool (*OrderedQ)(T a, T b)) {
	int p1,p2;
	p1=0;
	p2=0;
	while((p1<len1)&&(p2<len2)) {
		// cout << "p1,p2: " << p1 << " " << p2 << endl;
		if((*OrderedQ)(list1[p1],list2[p2])) {
			target[p1+p2]=list1[p1];
			p1++;
		} else {
			target[p1+p2]=list2[p2];
			p2++;
		}
	}
	while(p1<len1) {
		target[p1+p2]=list1[p1];
		p1++;
	}
	while(p2<len2) {
		target[p1+p2]=list2[p2];
		p2++;
	}

}


template<class T> void MSmergeSort(T* list1, int len, bool (*OrderedQ)(T a, T b)) {
	T *list2,*tList,*pList,*buffer;
	int res,i;
	list2=(T*) malloc(sizeof(T)*len);
	res=1;
	tList=list1;
	pList=list2;
	while(res<len) {
		//cout << "res: " << res << endl;
		// nBins=(int) ( (len-1)/res )+1 anzahl der bins mit länge maximal res
		// (int) ( (nBins-1)/2 )+1 anzahl der paare, aufgerundet, i.e. letztes paar ist evtl nur ein bin,
		//     aber MSmergeLists behandelt das korrekt
		for(i=0;i<(int)(((int) (len-1)/res)/2)+1;i++) {
			MSmergeLists<T>(tList+2*i*res, tList+(2*i+1)*res,  pList+2*i*res,
				min(res,len-(2*i)*res), // check if there is actually still res fields in that bin, otherwise use remaining fields
				max(0,min(res,len-(2*i+1)*res)), // same here, but here also catch, if the bin does not even exist any more
				OrderedQ);
		}
		res=res*2;
		buffer=tList;
		tList=pList;
		pList=buffer;
	}
	if(tList!=list1) {
		for(i=0;i<len;i++) {
			list1[i]=tList[i];
		}
	}
	free(list2);
}

