#include <bits/stdc++.h>
#pragma GCC target("avx2")
#include <immintrin.h>
using i32=int32_t;
using i64=int64_t;
using u32=uint32_t;
using u64=uint64_t;
#define IL __inline__ __attribute__((always_inline))
#define RC(T,x) reinterpret_cast<T>(x)
#define LG2(x) std::__lg(x)
#define CRZ(x) __builtin_ctzll(x)
namespace No_Poly{
using idt=std::size_t;
constexpr u32 M=998244353u;
constexpr idt pool_MB=64;
std::byte _mem_pool[pool_MB<<20],*_now=_mem_pool;
template<idt align>inline bool is_align(const void*mem){
	static_assert(std::has_single_bit(align));
	return (RC(idt,mem)&(align-1))==0;
}
template<idt align>inline void*to_align(void*mem){
	static_assert(std::has_single_bit(align));
	return RC(void*,idt((RC(idt,mem)+align-1)&(-align)));
}
template<class T>concept trivialT=std::is_trivial_v<T>;
struct pl_alcor{
    std::byte*t;
    pl_alcor():t(_now){}
    ~pl_alcor(){_now=t;}
    template<idt al>void*raw(idt n){
        void*p=to_align<al>(_now);
        return _now=((std::byte*)p)+n,p;
    }
};
template<trivialT T,idt al>struct wrap_pl_alcor:pl_alcor{
	T*operator ()(idt n){return (T*)pl_alcor::raw<al>(n*sizeof(T));}
};
using _alcr=wrap_pl_alcor<u32,32>;
template<class T,idt align>struct aligned_alcor{
	static_assert(std::has_single_bit(align));
	typedef T value_type;
	T*allocate(idt n){return new(std::align_val_t(align))T[n];}
	template<class U>struct rebind{using other=aligned_alcor<U,align>;};
	void deallocate(T*p,idt){::operator delete[](p,std::align_val_t(align));}
};
using vec=std::vector<u32,aligned_alcor<u32,32> >;

using u32x8=__attribute((vector_size(32))) u32;
using u64x4=__attribute((vector_size(32))) u64;
using I256=__m256i;
using I256u=__m256i_u;
template<bool align=true>IL u32x8 load(const void*data){
	if constexpr(align){return (u32x8)_mm256_load_si256((const I256*)data);}
	return (u32x8)_mm256_loadu_si256((const I256u*)data);
}
template<bool align=true>IL void store(const u32x8&x,void*data){
	if constexpr(align){return _mm256_store_si256((I256*)data,RC(I256,x));}
	return _mm256_storeu_si256((I256u*)data,RC(I256,x));
}
IL u64x4 fus_mul(const u32x8&x,const u32x8&y){return RC(u64x4,_mm256_mul_epu32(RC(I256,x),RC(I256,y)));}
IL u32x8 swaplohi128(const u32x8&x){return (u32x8)_mm256_permute2x128_si256(RC(I256,x),RC(I256,x),1);}
template<int typ>IL u32x8 shuffle(const u32x8&x){return RC(u32x8,_mm256_shuffle_epi32(RC(I256,x),typ));}
template<int typ>IL u32x8 blend(const u32x8&x,const u32x8&y){return RC(u32x8,_mm256_blend_epi32(RC(I256,x),RC(I256,y),typ));}
IL u32x8&x8(u32*data){return*((u32x8*)data);}
IL const u32x8&x8(const u32*data){return*((const u32x8*)data);}
IL u32x8 min_u32(const u32x8&x,const u32x8&y){return RC(u32x8,_mm256_min_epu32(RC(I256,x),RC(I256,y)));}
constexpr IL u32x8 padd(u32 x){return (u32x8){x,x,x,x,x,x,x,x};}

constexpr u32 get_nr(u32 M){u32 Iv=2u-M;for(int i=0;i<4;++i){Iv*=2-M*Iv;}return Iv;}
constexpr u32 pr_rt(u32 M){u32 qed=0,n=M-1,d[11]={};for(u32 i=2;i*i<=n;++i){if(n%i==0){d[qed++]=i;do{n/=i;}while(n%i==0);}}if(n>1){d[qed++]=n;}for(u32 g=2,r=0;;++g){for(u32 i=0;i<qed;++i){u32 b=(M-1)/d[i],a=g;for(r=1;b;b>>=1,a=u64(a)*a%M){b&1?r=u64(r)*a%M:r;}if(r==1){break;}}if(r!=1){return g;}}}
template<u32 M>constexpr u32 sf_md(i64 x){return x%=M,x<0?x+M:x;}
constexpr idt bcl(idt x){return ((x<2)?1:idt(2)<<LG2(x-1));}

constexpr u32 R=(-M)%M,E={},nR=M-R,M2=M*2,iv=get_nr(M),niv=-iv,R2=(-u64(M))%M;
constexpr u32 shrk(u32 x){return x<M?x:x-M;}
constexpr u32 dil2(u32 x){return x>>31?x+M2:x;}
constexpr u32 reduce(u64 x){return (x+u64(u32(x)*niv)*M)>>32;}
constexpr u32 reduce_s(u64 x){u32 r=(x>>32)-((u64(u32(x)*iv)*M)>>32);return r>>31?r+M:r;}

constexpr u32 add(u32 x,u32 y){return dil2(x+y-M2);}
constexpr u32 sub(u32 x,u32 y){return dil2(x-y);}
constexpr u32 mul(u32 x,u32 y){return reduce(u64(x)*y);}
constexpr u32 mul_s(u32 x,u32 y){return reduce_s(u64(x)*y);}
constexpr u32 qpw(u32 a,u32 b,u32 r=R){for(;b;b>>=1,a=mul(a,a)){b&1?r=mul(r,a):r;}return r;}
constexpr u32 inv(u32 x){return qpw(x,M-2);}
constexpr u32 dvs(u32 x,u32 y){return qpw(y,M-2,x);}
constexpr u32 neg(u32 x){return M2-x;}
constexpr u32 in(u32 x){return mul(x,R2);}
constexpr u32 in_s(u32 x){return mul_s(x,R2);}
constexpr u32 out(u32 x){u32 r=(x+(u64(x*niv)*M))>>32;return r<M?r:r-M;}
constexpr bool equals(u32 x,u32 y){return out(x)==out(y);}
constexpr void clr(u32&x){x=E;}

constexpr u32x8 Rx8=padd(R),Ex8=padd(E),Mx8=padd(M),M2x8=padd(M2),nivx8=padd(niv);
IL u32x8 shrk(const u32x8&x){return min_u32(x,x-M);}
IL u32x8 dil2(const u32x8&x){return min_u32(x,x+M2);}
IL u32x8 add(const u32x8&x,const u32x8&y){return dil2(x+y-M2x8);}
IL u32x8 sub(const u32x8&x,const u32x8&y){return dil2(x-y);}
IL u32x8 mul(const u32x8&x,const u32x8&y){
	u32x8 z=nivx8*x*y;
	return blend<0xaa>(RC(u32x8,(fus_mul(x,y)+fus_mul(z,Mx8))>>32),RC(u32x8,(fus_mul(u32x8(u64x4(x)>>32),u32x8(u64x4(y)>>32))+fus_mul(shuffle<0xf5>(z),Mx8))));
}
IL u32x8 qpw(const u32x8&y,u32 b,const u32x8&_r=Rx8){u32x8 x=y,r=_r;for(;b;x=mul(x,x),b>>=1){if(b&1){r=mul(r,x);}}return r;}
IL u32x8 inv(const u32x8&x){return qpw(x,M-2);}
IL u32x8 dvs(const u32x8&x,const u32x8&y){return qpw(y,M-2,x);}
IL u32x8 mul_s(const u32x8&x,const u32x8&y){return shrk(mul(x,y));}
IL u32x8 neg(const u32x8&x){return M2x8-x;}
IL void clr(u32x8&x){x=Ex8;}

constexpr u32 _Amul(u32 a,u32 b,u32 c){return mul(a+b,c);}
constexpr u32 _Smul(u32 a,u32 b,u32 c){return mul(a-b+M2,c);}
IL u32x8 _LMadd(const u32x8&x,const u32x8&y){return x+y;}
IL u32x8 _LMsub(const u32x8&x,const u32x8&y){return x-y+M2x8;}
IL u32x8 _LMnot(const u32x8&x){return min_u32(x,x-M2);}
IL u32x8 _Amul(const u32x8&a,const u32x8&b,const u32x8&c){return mul(a+b,c);}
IL u32x8 _Smul(const u32x8&a,const u32x8&b,const u32x8&c){return mul(a-b+M2x8,c);}
template<int typ>IL u32x8 Neg(const u32x8&x){return blend<typ>(x,M2x8-x);}
constexpr u32x8 powXx8(u32 X){
	u32 X2=mul_s(X,X),X3=mul_s(X2,X),X4=mul_s(X3,X),X5=mul_s(X4,X),X6=mul_s(X5,X),X7=mul_s(X6,X);
	return (u32x8){R,X,X2,X3,X4,X5,X6,X7};
}
constexpr u32 _ADmul(u32 a,u32 b,u32 c,u32 d){return reduce_s(u64(a)*b+u64(c)*d);}
IL u32x8 _ADmul(const u32x8&a,const u32x8&b,const u32x8&c,const u32x8&d){
	u32x8 z=nivx8*(a*b+c*d);
	return shrk(blend<0xaa>(RC(u32x8,(fus_mul(a,b)+fus_mul(c,d)+fus_mul(z,Mx8))>>32),RC(u32x8,(fus_mul(u32x8(u64x4(a)>>32),u32x8(u64x4(b)>>32))+fus_mul(u32x8(u64x4(c)>>32),u32x8(u64x4(d)>>32))+fus_mul(shuffle<0xf5>(z),Mx8)))));
}

constexpr u32 sqt(u32 x){
	u32 y=R,Im=E;
	while(shrk(qpw(Im=sub(mul(y,y),x),(M-1)>>1))==R){++y;}
	struct dZ{u32 r,i;};
	auto Mul=[Im](dZ a,dZ b)->dZ{return {_ADmul(a.r,b.r,mul(Im,a.i),b.i),_ADmul(a.r,b.i,a.i,b.r)};};
	dZ a={y,R},r={R,E};
	for(u32 b=(M+1)>>1;b;b>>=1,a=Mul(a,a)){if(b&1){r=Mul(r,a);}}
	u32 z=out(r.r);
	return in_s(std::min(z,M-z));
}

constexpr u32 Half=shrk(inv(in(2))),nHalf=M-Half;

struct vec_v{
    u32x8 v;
    constexpr vec_v(u32 x):v{padd(x)}{}
    constexpr operator u32()const{return v[0];}
    constexpr operator u32x8()const{return v;}
};

inline void vec_op(auto f,idt n,auto&&op){idt i=0;for(;i+7<n;i+=8){op(x8(f+i));}for(;i<n;++i){op(f[i]);}}
inline void vec_op(auto f,auto g,idt n,auto&&op){idt i=0;for(;i+7<n;i+=8){op(x8(f+i),x8(g+i));}for(;i<n;++i){op(f[i],g[i]);}}
inline void vec_op(auto f,auto g,auto h,idt n,auto&&op){idt i=0;for(;i+7<n;i+=8){op(x8(f+i),x8(g+i),x8(h+i));}for(;i<n;++i){op(f[i],g[i],h[i]);}}
inline void vec_op(auto f,auto g,auto h,auto o,idt n,auto&&op){idt i=0;for(;i+7<n;i+=8){op(x8(f+i),x8(g+i),x8(h+i),x8(o+i));}for(;i<n;++i){op(f[i],g[i],h[i],o[i]);}}

constexpr u32 _g=in_s(pr_rt(M));
namespace raw_ntt{
constexpr int lml=CRZ(M-1);
struct P_R_Tab{
	u32 t[lml+1];
	constexpr P_R_Tab(u32 G):t{}{
        t[lml]=shrk(qpw(G,(M-1)>>lml));
        for(int i=lml;i>0;--i){t[i-1]=mul_s(t[i],t[i]);}
    }
	constexpr u32 operator[](int i)const{return t[i];}
};
struct ntt_info_base4x8{
	u32 rt3[lml-2],rt3_I[lml-2];
	u32x8 rt4ix8[lml-3],rt4ix8_I[lml-3];
	constexpr ntt_info_base4x8(const P_R_Tab&w,const P_R_Tab&wI):rt3{},rt3_I{},rt4ix8{},rt4ix8_I{}{   
		u32 pr=R,pr_I=R;
		for(int i=0;i<lml-2;pr=mul(pr,wI[i+3]),pr_I=mul(pr_I,w[i+3]),++i){
			rt3[i]=mul_s(pr,w[i+3]),rt3_I[i]=mul_s(pr_I,wI[i+3]);
		}
		pr=R,pr_I=R;
		for(int i=0;i<lml-3;pr=mul(pr,wI[i+4]),pr_I=mul(pr_I,w[i+4]),++i){
			rt4ix8[i]=powXx8(mul_s(pr,w[i+4])),rt4ix8_I[i]=powXx8(mul_s(pr_I,wI[i+4]));
		}
	}
};
constexpr P_R_Tab rt1={_g},rt1_I={inv(_g)};
constexpr ntt_info_base4x8 iab4={rt1,rt1_I};
constexpr u32 Img=rt1[2];
constexpr u32x8 Imgx8=padd(Img);
template<bool strict=false>inline void dif_2(u32&x,u32&y){
	u32 sum=add(x,y),diff=sub(x,y);x=sum,y=diff;
	if constexpr(strict){x=shrk(x),y=shrk(y);}
}
template<bool strict=false>inline void dif_4(u32&x,u32&y,u32&z,u32&w){
	u32 a=sub(x,z),b=_Smul(y,w,Img);
	x=add(x,z),y=add(y,w),z=add(a,b),w=sub(a,b),a=add(x,y),b=sub(x,y),x=a,y=b;
	if constexpr(strict){x=shrk(x),y=shrk(y),z=shrk(z),w=shrk(w);}
}
template<bool strict=false>inline void vec_dif_base4(u32x8*f,idt n){
	idt L=n>>1;
	if(CRZ(n)&1){
		for(idt j=0;j<L;++j){
			auto x=f[j],y=f[j+L];
			f[j]=_LMadd(x,y),f[j+L]=_LMsub(x,y);
		}
		L>>=1;
	}
	L>>=1;
	for(idt l=L<<2,k;L;l=L,L>>=2){
		u32 r=R,r2=R,r3=nR;k=1;
		for(auto i=f;i!=(f+n);r=mul_s(r,iab4.rt3[CRZ(k++)]),r2=mul_s(r,r),r3=mul_s(r2,neg(r)),i+=l){
			auto rx8=padd(r),r2x8=padd(r2),r3x8=padd(r3);
			for(auto F0=i,F1=F0+L,F2=F1+L,F3=F2+L;F3!=i+l;++F0,++F1,++F2,++F3){
				auto f0=_LMnot(*F0),f1=mul(*F1,rx8),f2=mul(*F2,r2x8),f3=mul(*F3,r3x8);
				auto f1f3=_Amul(f1,f3,Imgx8),f02=add(f0,f2),f13=sub(f1,f3),f_02=sub(f0,f2);
				*F0=_LMadd(f02,f13),*F1=_LMsub(f02,f13),*F2=_LMadd(f_02,f1f3),*F3=_LMsub(f_02,f1f3);
			}
		}
	}
	constexpr u32x8 pr2={R,R,R,Img,R,R,R,Img},pr4={R,R,R,R,R,rt1[3],Img,mul_s(Img,rt1[3])};
	auto rx8=Rx8;
	for(idt i=0;i<n;rx8=mul_s(rx8,iab4.rt4ix8[CRZ(++i)])){
		auto&fi=f[i];fi=mul(fi,rx8);
		fi=_Amul(Neg<0xf0>(fi),swaplohi128(fi),pr4);
		fi=_Amul(Neg<0xcc>(fi),shuffle<0x4e>(fi),pr2);
		fi=sub(shuffle<0xb1>(fi),Neg<0x55>(fi));
		if constexpr(strict){fi=shrk(fi);}
	}
}
template<u32 fx>inline void dit_2(u32&x,u32&y){
	constexpr u32 iv2=mul_s(inv(in(2)),fx);
	u32 a=_Amul(x,y,iv2),b=_Smul(x,y,iv2);x=a,y=b;
}
template<u32 fx>inline void dit_4(u32&x,u32&y,u32&z,u32&w){
	constexpr u32 iv4=mul_s(inv(in(4)),fx),Imgi4=mul_s(iv4,Img);
	u32 a=_Amul(x,y,iv4),b=_Smul(x,y,iv4);
	x=a,y=b,a=_Amul(z,w,iv4),b=_Smul(w,z,Imgi4),z=sub(x,a),w=sub(y,b),x=add(x,a),y=add(y,b);
}
template<u32 fx>inline void vec_dit_base4(u32x8*f,idt n){
	idt L=1;
	constexpr u32 nR2=in_s(nR),M8=(M-1)>>3;
	constexpr u32x8 pr2={nR2,nR2,nR2,in(Img),nR2,nR2,nR2,in(Img)},pr4={fx,fx,fx,fx,fx,mul_s(fx,rt1_I[3]),mul_s(fx,rt1_I[2]),mul_s(fx,mul_s(rt1_I[2],rt1_I[3]))};
	auto rx8=padd(M8>>CRZ(n));
	for(idt i=0;i<n;rx8=mul_s(rx8,iab4.rt4ix8_I[CRZ(++i)])){
		auto&fi=f[i];
		fi=_Amul(Neg<0xaa>(fi),shuffle<0xb1>(fi),pr2);
		fi=_Amul(Neg<0xcc>(fi),shuffle<0x4e>(fi),pr4);
		fi=_Amul(Neg<0xf0>(fi),swaplohi128(fi),rx8);
	}
	for(idt l=L<<2,k;L<(n>>1);L=l,l<<=2){
		u32 r=R,r2=R,r3=R;k=1;
		for(auto i=f;i!=(f+n);r=mul_s(r,iab4.rt3_I[CRZ(k++)]),r2=mul_s(r,r),r3=mul_s(r2,r),i+=l){
			auto rx8=padd(r),r2x8=padd(r2),r3x8=padd(r3);
			for(auto F0=i,F1=F0+L,F2=F1+L,F3=F2+L;F3!=i+l;++F0,++F1,++F2,++F3){
				auto f0=*F0,f1=*F1,f2=neg(*F2),f3=*F3;
				auto f2f3=_Amul(f3,f2,Imgx8),f01=add(f0,f1),f23=sub(f2,f3),f_01=sub(f0,f1);
				*F0=sub(f01,f23),*F1=_Amul(f_01,f2f3,rx8),*F2=_Amul(f01,f23,r2x8),*F3=_Smul(f_01,f2f3,r3x8);
			}
		}
	}
	if(CRZ(n)&1){
		for(idt j=0;j<L;++j){
			auto x=f[j],y=f[j+L];
			f[j]=add(x,y),f[j+L]=sub(x,y);
		}
	}
}
template<bool strict=false>inline void dif(u32*A,idt lim){
	switch(lim){
		case 1:break;
		case 2:dif_2<strict>(A[0],A[1]);break;
		case 4:dif_4<strict>(A[0],A[1],A[2],A[3]);break;
		default:vec_dif_base4<strict>((u32x8*)A,lim>>3);
	}
}
template<u32 fx=R>inline void dit(u32*A,idt lim){
	switch(lim){
		case 1:if constexpr(!equals(fx,R)){A[0]=mul(A[0],fx);}break;
		case 2:dit_2<fx>(A[0],A[1]);break;
		case 4:dit_4<fx>(A[0],A[1],A[2],A[3]);break;
		default:vec_dit_base4<fx>((u32x8*)A,lim>>3);
	}
}
}//namespace raw_ntt
using raw_ntt::dif;
using raw_ntt::dit;
//u32* u32 const u32* idt 
template<trivialT T>inline T*cpy(T*f,const T*g,idt n){return (T*)memcpy(f,g,n*sizeof(T));}
template<trivialT T>inline T*clr(T*f,idt n){return (T*)memset(f,0,n*sizeof(T));}
template<trivialT T>inline T*rcpy(T*f,const T*g,idt n){return std::reverse_copy(g,g+n,f),f;}
void dot(u32*f,const u32*g,idt n){vec_op(f,g,n,[](auto&fi,auto&gi){fi=mul(fi,gi);});}
void dot(u32*f,const u32*g,const u32*h,idt n){vec_op(f,g,h,n,[](auto&fi,auto&gi,auto&hi){return fi=mul(gi,hi);});}
void add(u32*f,const u32*g,idt n){vec_op(f,g,n,[](auto&fi,auto&gi){fi=add(fi,gi);});}
void add(u32*f,const u32*g,const u32*h,idt n){vec_op(f,g,h,n,[](auto&fi,auto&gi,auto&hi){return fi=add(gi,hi);});}
void sub(u32*f,const u32*g,idt n){vec_op(f,g,n,[](auto&fi,auto&gi){fi=sub(fi,gi);});}
void sub(u32*f,const u32*g,const u32*h,idt n){vec_op(f,g,h,n,[](auto&fi,auto&gi,auto&hi){return fi=sub(gi,hi);});}
void adddot(u32*a,const u32*b,const u32*c,const u32*d,idt n){vec_op(a,b,c,d,n,[](auto&A,auto&B,auto&C,auto&D){A=_ADmul(A,B,C,D);});}
void vec_multi_iv(u32x8*f,const u32x8*g,idt n){
	if(n==0){return;}
	f[0]=g[0];for(idt i=1;i<n;++i){f[i]=mul(f[i-1],g[i]);}
	f[n-1]=inv(f[n-1]);
	for(idt i=n-1;i;--i){
		u32x8 ivi=f[i];
		f[i]=mul(ivi,f[i-1]),f[i-1]=mul(ivi,g[i]);
	}
}
void multi_iv(u32*f,const u32*g,idt n){
	vec_multi_iv((u32x8*)f,(u32x8*)g,n>>3);
	for(idt i=(n>>3)<<3;i<n;++i){f[i]=inv(g[i]);}
}
template<u32 fx=R>void conv(u32*f,u32*g,idt lim){dif(f,lim),dif(g,lim),dot(f,g,lim),dit<fx>(f,lim);}

template<class T>concept like_vec=requires(T a,idt b){a.data();a.size();a[b];a.resize(b);};
template<like_vec V,idt alz=16>struct sim_inf_seq{
	static_assert(std::has_single_bit(alz));
	using T=V::value_type;
	std::function<void(T*,idt,idt)> f;
	mutable V v;
	sim_inf_seq(auto&&F):f{F}{}
	T*rsv(idt l)const{
		idt ol=v.size();
		if(l>ol)[[unlikely]]{l=std::max((l+alz-1)&-alz,ol<<1),v.resize(l),f(v.data(),ol,l);}
		return v.data();
	}
	const T&operator[](idt pos)const{return rsv(pos+1)[pos];}
};
sim_inf_seq<vec> 
Id=[](u32*f,idt l,idt r){
	for(;l<8;++l){f[l]=in(l);}
	constexpr auto va8=padd(in(8));
	for(auto i=l;i<r;i+=8){x8(f+i)=add(x8(f+i-8),va8);}
},
Iv=[](u32*f,idt l,idt r){
	auto id=Id.rsv(r);
	if(l<8){x8(f)=inv(x8(id)),l=8;}
	vec_multi_iv((u32x8*)(f+l),(const u32x8*)(id+l),(r-l)>>3);
},
Fac=[](u32*f,idt l,idt r){
	auto id=Id.rsv(r);
	if(l==0){f[0]=R,l=1;}
	for(auto i=l;i<r;++i){f[i]=mul(f[i-1],id[i]);}
},
IFac=[](u32*f,idt l,idt r){
	auto iv=Iv.rsv(r);
	if(l==0){f[0]=R,l=1;}
	for(auto i=l;i<r;++i){f[i]=mul(f[i-1],iv[i]);}
};

void deriv(u32*f,const u32*g,idt n){
	idt i=1;
	auto id=Id.rsv(n);
	if(n>16){
		for(;i<8;++i){f[i-1]=mul(g[i],id[i]);}
		for(;i+7<n;i+=8){store<false>(mul(load(g+i),x8(id+i)),f+i-1);}
	}
	for(;i<n;++i){f[i-1]=mul(g[i],id[i]);}
}
void integ(u32*f,const u32*g,idt n,u32 C=E){
	idt i=n;
	auto iv=Iv.rsv(n);
	if(n>16){
		for(;(i&7)!=7;--i){f[i]=mul(g[i-1],iv[i]);}
		for(i-=7;i>0;i-=8){store(mul(load<false>(g+i-1),x8(iv+i)),f+i);}i+=7;
	}
	for(;i>0;--i){f[i]=mul(g[i-1],iv[i]);}f[0]=C;
}

void scan(u32*f,idt n=1){
	for(idt i=0;i<n;++i){
		std::cin>>f[i],f[i]=in(f[i]);
	}
}
void print(const u32*f,idt n=1){
	for(idt i=0;i<n;++i){
		std::cout<<out(f[i])<<" \n"[i+1==n];
	}
}

void inv(u32*f,const u32*g,idt n){
	auto lim=bcl(n);
	_alcr alc;
	auto o=alc(lim),h=alc(lim);
	f[0]=inv(g[0]);
	for(idt t=2,m=1,xl;t<=lim;m=t,t<<=1){
		xl=std::min(n,t),clr(cpy(o,g,xl)+xl,t-xl),clr(cpy(h,f,m)+m,m),conv(o,h,t);
		clr(o,m),dif(o,t),dot(o,h,t),dit<nR>(o,t),cpy(f+m,o+m,xl-m);
	}
}
void quo(u32*f,const u32*g,const u32*h,idt n){
	_alcr alc;
	if(n<=64){
		idt lim=bcl(n<<1);auto o=alc(lim),s=alc(lim);
		inv(o,h,n),clr(o+n,lim-n),cpy(s,g,n),clr(s+n,lim-n),conv(o,s,lim),cpy(f,o,n);return;
	}
	idt bn=bcl(n)>>4,bt=(n+bn-1)/bn,bn2=bn<<1;
	auto o=alc(bn2),A=alc(bn2);
	inv(o,h,bn),clr(o+bn,bn),clr(cpy(A,g,bn)+bn,bn),conv(A,o,bn2);
	auto nh=alc(bn2*bt),Nf=alc(bn2*(bt-1)),nf=Nf;
	cpy(f,A,bn),clr(cpy(nh,h,bn)+bn,bn),dif(nh,bn2);
	for(idt ds=bn,xl;ds<n;ds+=bn){
		xl=std::min(bn,n-ds),nh+=bn2;clr(cpy(nh,h+ds,xl)+xl,bn2-xl),dif(nh,bn2);
		clr(cpy(nf,f+ds-bn,bn)+bn,bn),dif<1>(nf,bn2),clr(A,bn2),nf+=bn2;auto nH=nh,nF=Nf,nH1=nH-bn2;
		for(idt dj=0;dj<ds;dj+=bn,nH-=bn2,nH1-=bn2,nF+=bn2){
			for(idt i=0;i<bn;i+=8){x8(A+i)=sub(x8(A+i),_Amul(x8(nH+i),x8(nH1+i),x8(nF+i)));}
			for(idt i=bn;i<bn2;i+=8){x8(A+i)=sub(x8(A+i),_Smul(x8(nH+i),x8(nH1+i),x8(nF+i)));}
		}
		dit(A,bn2),clr(A+bn,bn),add(A,g+ds,xl),dif(A,bn2),dot(A,o,bn2),dit(A,bn2),cpy(f+ds,A,xl);
	}
}
void dvs(u32*q,const u32*f,const u32*g,idt n,idt m){
	idt lm=n-m+1,R=std::min(m,lm);
	_alcr alc;
	auto o=alc(lm);
	clr(rcpy(o,g+m-R,R)+R,lm-R),quo(q,rcpy(q,f+m-1,lm),o,lm),std::reverse(q,q+lm);
}
void dvs(u32*q,u32*r,const u32*f,const u32*g,idt n,idt m){
	dvs(q,f,g,n,m);
	idt lm=bcl(std::min(n,m+m-3)),u=m-1,v=std::min(u,n-u);
	if(v<=16){
		for(idt i=0,k;i<u;++i){
			for(k=0,r[i]=f[i];k<std::min(v,i+1);++k){r[i]=sub(r[i],mul(q[k],g[i-k]));}
		}
	}
	else{
		_alcr alc;
		auto o=alc(lm),h=alc(lm);
		clr(cpy(o,g,u)+u,lm-u),clr(cpy(h,q,v)+v,lm-v),conv(o,h,lm),sub(r,f,o,u);
	}
}
void ln(u32*f,const u32*g,idt n){dot(f,Id.rsv(n),g,n),quo(f,f,g,n),dot(f,Iv.rsv(n),n);}
template<bool c_inv>void __expi(u32*f,u32*h,const u32*g,idt n){
	f[0]=h[0]=R;if(n==1){return;}
	_alcr alc;
	auto lim=bcl(n);
	auto id=Id.rsv(lim),iv=Iv.rsv(lim);
	auto o=alc(lim),A=alc(lim),B=alc(lim);
	clr(A,lim),A[0]=A[1]=R;
	for(idt t=2,m=1,xl;t<=lim;m=t,t<<=1){
		xl=std::min(n,t),dot(o,id,g,m),dif(o,m),dot(o,A,m),dit(o,m);dot(o+m,f,id,m);
		vec_op(o+m,o,m,[](auto&fi,auto&gi){fi=sub(fi,gi),clr(gi);}),dif(o,t);
		clr(cpy(B,h,m)+m,m),dif(B,t),dot(o,B,t),dit(o,t),dot(clr(o,m)+m,iv+m,m);
		sub(o+m,g+m,xl-m),dif(o,t),dot(A,o,t),dit<nR>(A,t),cpy(f+m,A+m,xl-m);
		if(c_inv||(t!=lim)){
			cpy(A,f,m),dif(A,std::min(t<<1,lim)),dot(o,A,B,t),dit(o,t),clr(o,m);
			dif(o,t),dot(o,B,t),dit<nR>(o,t),cpy(h+m,o+m,xl-m);
		}
	}
}
void exp(u32*f,const u32*g,idt n){
	_alcr alc;
	if(n<=64){return __expi<false>(f,alc(n),g,n);}
	idt bn=bcl(n)>>4,bt=(n+bn-1)/bn,bn2=bn<<1;
	auto o=alc(bn2),h=alc(bn2);
	auto id=Id.rsv(n),iv=Iv.rsv(n);
	__expi<true>(f,h,g,bn),clr(h+bn,bn),dif(h,bn2);
	auto ng=alc(bn2*bt),Nf=alc(bn2*(bt-1)),nf=Nf;
	dot(ng,g,id,bn),clr(ng+bn,bn),dif(ng,bn2);
	for(idt ds=bn,xl;ds<n;ds+=bn){
		xl=std::min(bn,n-ds),ng+=bn2,dot(ng,g+ds,id+ds,xl),clr(ng+xl,bn2-xl),dif(ng,bn2);
		clr(cpy(nf,f+ds-bn,bn)+bn,bn),dif<1>(nf,bn2),clr(o,bn2),nf+=bn2;
		auto nG=ng,nF=Nf,nG1=nG-bn2;
		for(idt dj=0;dj<ds;dj+=bn,nG-=bn2,nG1-=bn2,nF+=bn2){
			for(idt i=0;i<bn;i+=8){x8(o+i)=sub(x8(o+i),_Amul(x8(nG+i),x8(nG1+i),x8(nF+i)));}
			for(idt i=bn;i<bn2;i+=8){x8(o+i)=sub(x8(o+i),_Smul(x8(nG+i),x8(nG1+i),x8(nF+i)));}
		}
		dit(o,bn2),clr(o+bn,bn),dif(o,bn2),dot(o,h,bn2),dit<nR>(o,bn2);
		dot(o,iv+ds,xl),clr(o+xl,bn2-xl),dif(o,bn2),dot(o,Nf,bn2),dit(o,bn2),cpy(f+ds,o,xl);
	}
}
void sqrtinv(u32*f,const u32*g,idt n){
	auto lim=bcl(n);
	_alcr alc;
	auto o=alc(lim*2),h=alc(lim*2);
	f[0]=inv(sqt(g[0]));
	for(idt r=4,t=2,m=1,xl;t<=lim;m=t,t=r,r<<=1){
		xl=std::min(n,t),clr(cpy(o,f,m)+m,r-m),clr(cpy(h,g,xl)+xl,r-xl),dif(o,r),dif(h,r);
		vec_op(h,o,r,[](auto&fi,auto&gi){fi=mul(mul(fi,gi),mul(gi,gi));}),dit<nHalf>(h,r),cpy(f+m,h+m,xl-m);
	}
}
template<bool c_inv>void __sqrti(u32*f,u32*h,const u32*g,idt n){
	auto lim=bcl(n);
	f[0]=sqt(g[0]),h[0]=inv(f[0]);
	_alcr alc;
	auto o=alc(lim),H=alc(lim),F=alc(lim);
	F[0]=f[0];
	for(idt t=2,m=1,xl;t<=lim;m=t,t<<=1){
		xl=std::min(t,n),dot(F,F,m),dit(F,m);
		vec_op(F,F+m,g,g+m,m,[](auto&a0,auto&a1,auto&b0,auto&b1){a1=sub(sub(a0,b0),b1),clr(a0);});
		clr(cpy(H,h,m)+m,m),conv<nHalf>(F,H,t),cpy(f+m,F+m,xl-m);
		if(c_inv||(t!=lim)){
			dif(cpy(o,f,t),t),cpy(F,o,t),dot(o,H,t),dit(o,t),dif(clr(o,m),t),dot(o,H,t),dit<nR>(o,t),cpy(h+m,o+m,xl-m);
		}
	}
}
void sqrt(u32*f,const u32*g,idt n){
	_alcr alc;
	if(n<=64){return __sqrti<false>(f,alc(n),g,n);}
	idt bn=bcl(n)>>4,bt=(n+bn-1)/bn,bn2=bn<<1;
	auto o=alc(bn2),jok=alc(bn2);
	__sqrti<true>(f,o,g,bn),clr(o+bn,bn),dif(o,bn2);
	auto nf=alc(bn2*(bt-1)),Nf=nf;
	for(idt ds=bn,xl;ds<n;ds+=bn){
		xl=std::min(bn,n-ds),clr(cpy(nf,f+ds-bn,bn)+bn,bn),dif<1>(nf,bn2),nf+=bn2;
		auto nF=nf,nF1=nf-bn2,NF=Nf;
		for(idt i=0;i<bn;i+=8){x8(jok+i)=neg(mul(x8(nF1+i),x8(NF+i)));}
		for(idt i=bn;i<bn2;i+=8){x8(jok+i)=mul(x8(nF1+i),x8(NF+i));}
		for(idt dj=bn;nF-=bn2,nF1-=bn2,NF+=bn2,dj<ds;dj+=bn){
			for(idt i=0;i<bn;i+=8){x8(jok+i)=sub(x8(jok+i),_Amul(x8(nF1+i),x8(nF+i),x8(NF+i)));}
			for(idt i=bn;i<bn2;i+=8){x8(jok+i)=sub(x8(jok+i),_Smul(x8(nF+i),x8(nF1+i),x8(NF+i)));}
		}
		dit<nR>(jok,bn2),clr(jok+bn,bn),sub(jok,g+ds,xl),dif(jok,bn2),dot(jok,o,bn2),dit<nHalf>(jok,bn2),cpy(f+ds,jok,xl);
	}
}

void Ci(u32*f,u32 z,idt n){
	u32 x=R;idt i=0;
	for(;i<std::min<idt>(n,8);++i,x=mul(x,z)){f[i]=x;}
	if(n>16){
		u32x8 xx8=padd(x);
		for(;i+7<n;i+=8){x8(f+i)=mul(x8(f+i-8),xx8);}
	}
	for(;i<n;++i){f[i]=mul(f[i-1],z);}
}
void mulk(u32 k,u32*f,const u32*g,idt n){
	auto k_v=vec_v{k};
	vec_op(f,g,n,[k_v](auto&fi,auto&gi){fi=mul(gi,k_v);});
}
namespace ntt_op{
using namespace raw_ntt;
sim_inf_seq<vec> 
sim_wn=[](u32*f,idt l,idt r){
	for(idt i=std::max<idt>(1,l);i<r;i<<=1){
		Ci(f+i,rt1[LG2(i)+1],i);
	}
},
bv_wn=[](u32*f,idt l,idt r){
	if(l==0){f[0]=R,l=1;}
	for(idt i=l;i<r;i<<=1){mulk(rt1[LG2(i)+1],f+i,f,i);}
};
template<bool strict=false>void dif2(u32*f,idt l){
	cpy(f+l,f,l),dit(f+l,l),dot(f+l,sim_wn.v.data()+l,l),dif<strict>(f+l,l);
}
constexpr u32 Two=in(2);
template<bool strict=false>void dif2_c1(u32*f,idt l){
	cpy(f+l,f,l),dit(f+l,l),dot(f+l,sim_wn.v.data()+l,l),f[l]=sub(Two,f[l]),dif<strict>(f+l,l);
}
template<bool strict=false>void dif2_xn(u32*f,idt l){
	cpy(f+l,f,l),dit(f+l,l),dot(f+l,sim_wn.v.data()+l,l),f[l]=sub(f[l],Two),dif<strict>(f+l,l);
}
void rev_dif(u32*f,const u32*g,idt l){
	dot(f,bv_wn.v.data(),g,l);
	for(idt i=2;i<l;i<<=1){std::reverse(f+i,f+i+i);}
}
idt locate_wn(u32 w){
	if(shrk(w)==R){return 0;}
	idt res=locate_wn(mul(w,w))<<1;
	return res+(shrk(w)!=shrk(bv_wn.v[res]));
}
}//namespace ntt_op
u32 eval(u32 x,const u32*f,idt n){
	u32 xn=R,res=E;
	for(idt i=0;i<n;++i){res=add(res,mul(xn,f[i])),xn=mul(xn,x);}
	return res;
}
void eval(u32*res,const u32*f,const u32*o,idt n,idt m){
	if(std::min(n,m)<=16){
		for(idt i=0;i<m;++i){res[i]=eval(o[i],f,n);}
		return;
	}
	using namespace ntt_op;
	static u32*GG[lml];
	_alcr alc;
	idt lm=bcl(std::max(n,m)),lm2=lm*2,m2=0;
	int lgn=std::__lg(lm);
	auto buf=alc(lm*3),pwh=alc(m),nw=GG[0]=alc(m*2),lt=nw;
	clr(cpy(buf,f,n)+n,lm-n),dif(buf,lm),bv_wn.rsv(lm);
	vec_op(pwh,o,m,[lgn](auto&fi,auto&gi){
        constexpr vec_v RR={R};
		fi=gi;
		for(int i=0;i<lgn;++i){fi=mul(fi,fi);}
		fi=shrk(sub(RR,fi));
	});
	for(idt i=0;i<m;++i){
		if(pwh[i]){nw[m2]=sub(R,o[i]),nw[m2|1]=sub(o[i],nR),m2+=2;}
        else{res[i]=buf[locate_wn(o[i])];}
	}
	if(m2>32){
		rev_dif(buf+lm,buf,lm),sim_wn.rsv(lm);
		for(int dep=1;dep<lgn;++dep,lt=nw){
			idt t=idt(1)<<dep,t2=t<<1,l=0,r=t;
			nw=GG[dep]=alc((m2+t2-1)&-t2);
			for(;r<m2;l+=t2,r+=t2){dot(nw+l,lt+l,lt+r,t),dif2_c1(nw+l,t);}
			if(l<m2){cpy(nw+l,lt+l,t),dif2(nw+l,t);}
		}
		dot(buf,lt,lt+lm,lm),multi_iv(buf+lm2,buf,lm),dot(buf+lm,buf+lm2,lm);
		dot(buf,buf+lm,lt+lm,lm),dot(buf+lm,lt,lm),dit(buf,lm),dit(buf+lm,lm);
		for(int dep=lgn-1;lt=GG[dep-1],dep>0;--dep){
			idt t=idt(1)<<dep,t2=t<<1,l=0,r=t,mid=t>>1;
			for(;r<m2;l+=t2,r+=t2){dif(buf+r,t),dot(buf+l,buf+r,lt+r,t),dot(buf+r,lt+l,t),dit(buf+l,t),dit(buf+r,t);}
			if(l<m2){cpy(buf+l+mid,buf+r+mid,mid);}
		}
		for(idt i=0,j=1;i<m;++i){if(pwh[i]){res[i]=mul(pwh[i],buf[j]),j+=2;}}
	}
	else{
		for(idt i=0;i<m;++i){if(pwh[i]){res[i]=eval(o[i],f,n);}}
	}
}
void intpol(u32*f,const u32*x,const u32*y,idt n){
	if(n==1){
		*f=*y;
		return;
	}
	using namespace ntt_op;
	static u32*GG[lml];
	_alcr alc;
	idt lm=bcl(n),lm2=lm*2,n2=n*2;
	int lgn=LG2(lm);
	sim_wn.rsv(lm);
	auto nw=GG[0]=alc(n2),lt=nw,buf=alc(lm*3);
	for(idt i=0;i<n;++i){nw[i<<1]=sub(R,x[i]),nw[i<<1|1]=sub(nR,x[i]);}
	for(int dep=1;dep<lgn;++dep,lt=nw){
		idt t=idt(1)<<dep,t2=t<<1,ed=n2&-t2,i=0;
		nw=GG[dep]=alc((n2+t2-1)&-t2);
		for(;i<ed;i+=t2){dot(nw+i,lt+i,lt+i+t,t),dif2_xn(nw+i,t);}
		if(i<n2){((n2-i)>t)?dot(nw+i,lt+i,lt+i+t,t):((void)cpy(nw+i,lt+i,t));dif2(nw+i,t);}
	}
	dot(buf,lt,lt+lm,lm),dit(buf,lm);
	if(n==lm){buf[lm]=R,buf[0]=sub(buf[0],R);}
	deriv(buf+lm2,buf,n+1),std::reverse(buf,buf+n+1),std::reverse(buf+lm2,buf+lm2+n);
	quo(buf+lm2,buf+lm2,buf,n),clr(buf+lm,lm-n),std::reverse_copy(buf+lm2,buf+lm2+n,buf+lm2-n);
	for(int dep=lgn;dep>0;--dep){
		lt=GG[dep-1];
		idt t=idt(1)<<dep,t2=t<<1,l=0,r=t,mid=t>>1;
		for(;r<n2;l+=t2,r+=t2){dif(buf+r,t),dot(buf+l,buf+r,lt+r,t),dit(buf+l,t),dot(buf+r,lt+l,t),dit(buf+r,t);}
		if(l<n2){cpy(buf+l+mid,buf+r+mid,mid);}
	}
	for(idt i=0;i<n;++i){buf[i]=buf[i<<1|1];}
	multi_iv(buf+lm2,buf,n);
	for(idt i=0;i<n;++i){buf[i<<1]=buf[i<<1|1]=mul_s(buf[i|lm2],y[i]);}
	for(int dep=1;lt=GG[dep-1],dep<lgn;++dep){
		idt t=idt(1)<<dep,t2=t<<1,l=0,r=t;
		for(;r<n2;l+=t2,r+=t2){adddot(buf+l,lt+r,buf+r,lt+l,t),dif2<true>(buf+l,t);}
		if(l<n2){dif2<true>(buf+l,t);}
	}
	adddot(buf,lt+lm,buf+lm,lt,lm),dit(buf,lm),cpy(f,buf,n);
}
void Neg(u32*f,const u32*g,idt l,idt r){
	for(;(l&7)&&(l<r);++l){f[l]=neg(g[l]);}
	vec_op(f+l,g+l,r-l,[](auto&fi,auto&gi){fi=neg(gi);});
}

template<class T>concept like_span=requires(T a,u32 b){b=*a.begin();a.end();a.size();};
template<class T>concept _mkseq_func=requires(T a,idt i,u32 b){b=a(i);};
template<class T>class poly_t:private T{
	private:
	u32*p;
	idt sz,cp;
	struct uinit0{idt x;};
	poly_t radd(const poly_t&b)const{
		poly_t r{uinit0{sz}};
		add(r.p,p,b.p,b.sz),cpy(r.p,p+b.sz,sz-b.sz),r.shrk();
		return r;
	}
    void shrk(){
        while(sz>0&&equals(p[sz-1],E)){--sz;}
    }
	explicit poly_t(uinit0 x):p{T::allocate(x.x)},sz{x.x},cp{x.x}{
        //std::clog<<"uinit0 construct\n";
    }
	explicit poly_t(idt n):p{T::allocate(n)},sz{n},cp{n}{
        //std::clog<<"init0 construct\n";
        clr(p,n);
    }
	poly_t rresize(idt n)const{
		poly_t r{uinit0{n}};
		clr(cpy(r.p,p,sz)+sz,n-sz);
		return r;
	}
	poly_t rmul(const poly_t&b)const{
		if(b.empty()){return {};}
		idt u=sz+b.sz-1;
		if(b.sz<=16){
			if(b.sz==1){return (*this)*b[0];}
			poly_t r(u);
			for(idt i=0;i<sz;++i){
				for(idt j=0;j<b.sz;++j){r[i+j]=add(r[i+j],mul(p[i],b[j]));}
			}
			return r;
		}
		idt lm=bcl(u);
		poly_t r=rresize(lm);
		_alcr alc;
		auto f=alc(lm);
		clr(cpy(f,b.p,b.sz)+b.sz,lm-b.sz),conv(r.p,f,lm),r.sz=u;
		return r;
	}
	void resize(idt nsz){
		if(nsz>sz){reserve(nsz),clr(p+sz,nsz-sz);}
		sz=nsz;
	}
	public:
	poly_t():p{nullptr},sz{0},cp{0}{}
	~poly_t(){T::deallocate(p,sz);}
	void clear(){sz=0;}
	poly_t(std::initializer_list<u32> lis):p(T::allocate(lis.size())),sz{lis.size()},cp{sz}{
		//std::clog<<"initializer_list construct\n";
		std::copy(lis.begin(),lis.end(),p);
	}
	void swap(poly_t&x){std::swap(p,x.p),std::swap(sz,x.sz),std::swap(cp,x.cp);}
	poly_t(const poly_t&x):p{T::allocate(x.sz)},sz{x.sz},cp{sz}{
        //std::clog<<"copy construct\n";
        cpy(p,x.p,sz);
    }
	poly_t(poly_t&&x):p{nullptr},sz{0},cp{0}{
        //std::clog<<"move construct\n";
        this->swap(x);
    }
	poly_t&operator=(poly_t&&x){
        //std::clog<<"move assign\n";
        this->swap(x);return*this;
    }
	poly_t&operator=(const poly_t&x){
        //std::clog<<"copy assign\n";
		if(p==x.p){return*this;}
		if(x.sz>cp){T::deallocate(p,sz),p=T::allocate(x.sz),cp=x.sz;}
		cpy(p,x.p,sz=x.sz);
		return*this;
	}
	template<like_span S>poly_t(const S&s):p{T::allocate(s.size())},sz{s.size()},cp{sz}{
		std::copy(s.begin(),s.end(),p),shrk();
	}
	u32&operator[](idt pos){return p[pos];}
	u32 operator[](idt pos)const{return p[pos];}
	u32 at(idt pos)const{return pos>=sz?0:p[pos];}
	void reserve(idt ncp){
		if(ncp>cp){
			u32*q=T::allocate(ncp);
			cpy(q,p,sz),T::deallocate(p,cp),cp=ncp,p=q;
		}
	}
	idt size()const{return sz;}
	idt capacity()const{return cp;}
	bool empty()const{return sz==0;}
	idt deg()const{return sz-1;}
	poly_t operator-()const{
		poly_t r{uinit0{sz}};
		vec_op(r.p,p,sz,[](auto&fi,auto&gi){fi=neg(gi);});
		return r;
	}
	poly_t&operator+=(const poly_t&b){
		if(b.sz>sz){reserve(b.sz),add(p,b.p,sz),cpy(p+sz,b.p+sz,b.sz-sz),shrk();}
		else{add(p,b.p,b.sz);}return*this;
	}
	friend poly_t operator+(const poly_t&a,const poly_t&b){
		return a.sz>b.sz?a.radd(b):b.radd(a);
	}
	poly_t&operator-=(const poly_t&b){
		if(b.sz>sz){reserve(b.sz),sub(p,b.p,sz),Neg(p,b.p,sz,b.sz),shrk();}
		else{sub(p,b.p,b.sz);}return*this;
	}
	friend poly_t operator-(const poly_t&a,const poly_t&b){
		poly_t r{uinit0{std::max(a.sz,b.sz)}};
		if(a.sz>b.sz){
			sub(r.p,a.p,b.p,b.sz),cpy(r.p+b.sz,a.p+b.sz,a.sz-b.sz);
		}
		else{
			sub(r.p,a.p,b.p,a.sz),Neg(r.p,b.p,a.sz,b.sz),r.shrk();
		}
		return r;
	}
	poly_t rsquare()const{
		idt u=sz+sz-1,lm=bcl(u);
		poly_t r=rresize(lm);
		dif(r.p,lm),dot(r.p,r.p,lm),dit(r.p,lm),r.sz=u;
		return r;
	}
	void square(){
		idt u=sz+sz-1,lm=bcl(u);
		resize(lm),dif(p,lm),dot(p,p,lm),dit(p,lm),sz=u;
	}
	poly_t&operator*=(const poly_t&b){
		auto t=std::min(sz,b.sz);
		if(t==0){return sz=0,*this;}
		if(p==b.p){
			square();
			return*this;
		}
		if(t<=16){
			idt u=sz+b.sz-1,z=sz;
			resize(u);
			for(idt i=z-1;~i;--i){
				u32 x=p[i];p[i]=E;
				for(idt j=0;j<b.sz;++j){p[i+j]=add(p[i+j],mul(x,b[j]));}
			}
			return*this;
		}
		idt u=sz+b.sz-1,lm=bcl(u);
		_alcr alc;
		auto f=alc(lm);
		cpy(f,b.p,b.sz),clr(f+b.sz,lm-b.sz),resize(lm),conv(p,f,lm),sz=u;
		return*this;
	}
	friend poly_t operator*(const poly_t&a,const poly_t&b){
		return a.sz<b.sz?b.rmul(a):a.rmul(b);
	}
	poly_t&operator*=(u32 k){
		mulk(k,p,p,sz);
		return*this;
	}
	friend poly_t operator*(const poly_t&f,u32 k){
		poly_t r{uinit0{f.sz}};
		mulk(k,r.p,f.p,f.sz);
		return r;
	}
	friend poly_t operator*(u32 k,const poly_t&f){
		return f*k;
	}
    std::pair<poly_t,poly_t> div_mod(const poly_t&b)const{
        assert(!b.empty());
        if(sz<b.sz){return {{},*this};}
        if(b.sz==1){return {(*this)*b[0],{}};}
        poly_t q{uinit0{sz-b.sz+1}},r{uinit0{b.sz-1}};
        dvs(q.p,r.p,p,b.p,sz,b.sz),r.shrk();
        return {std::move(q),std::move(r)};
    }
    friend poly_t operator/(const poly_t&a,const poly_t&b){
        assert(!b.empty());
        if(a.sz<b.sz){return {};}
        if(b.sz==1){return a*b[0];}
        poly_t q{uinit0{a.sz-b.sz+1}};
        dvs(q.p,a.p,b.p,a.sz,b.sz);
        return q;
    }
	poly_t&operator/=(const poly_t&b){
		return (*this)=(*this)/b;
	}
    friend poly_t operator%(const poly_t&a,const poly_t&b){
        return a.div_mod(b).second;
    }
	poly_t&operator%=(const poly_t&b){
		return (*this)=(*this)%b;
	}
    friend bool operator==(const poly_t&a,const poly_t&b){
        if(a.sz!=b.sz){return false;}
        for(idt i=0;i<a.sz;++i){
            if(!equals(a[i],b[i])){return false;}
        }
        return true;
    }
    friend bool operator!=(const poly_t&a,const poly_t&b){
        return !(a==b);
    }
    template<class _is=std::istream>void read(idt n,_is&is=std::cin){
        reserve(n);
        for(idt i=0;i<n;++i){is>>p[i],p[i]=in(p[i]);}
        sz=n,shrk();
    }
    template<class _os=std::ostream>void write(_os&os=std::cout)const{
        for(idt i=0;i<sz;++i){
            os<<out(p[i])<<" \n"[i+1==sz];
        }
    }
	template<class _os=std::ostream>void write(idt n,_os&os=std::cout)const{
        for(idt i=0;i<n;++i){
            os<<out(at(i))<<" \n"[i+1==n];
        }
    }
	void trunc(idt n){
		sz=std::min(sz,n);
	}
	poly_t rtrunc(idt n)const{
		return make_s(p,p+std::min(n,sz));
	}
	void pow(u32 n){
		if(n==0){*this={R};return;}
		if(n==1){return;}
		if(n==2){return square();}
		idt u=n*sz-n+1,lm=bcl(u);
		resize(lm),dif(p,lm),vec_op(p,[n](auto&fi){fi=qpw(fi,n);}),dit(p,lm),sz=u;
	}
	poly_t rpow(u32 n)const{
		if(n==0){return {R};}
		if(n==1){return *this;}
		if(n==2){return rsquare();}
		idt u=n*sz-n+1,lm=bcl(u);
		poly_t r=rresize(lm);
		dif(r.p,lm),vec_op(r.p,[n](auto&fi){fi=qpw(fi,n);}),dit(r.p,lm),r.sz=u;
		return r;
	}
	template<_mkseq_func F>static poly_t make_f(idt n,F&&op){
		poly_t r{uinit0{n}};
		for(idt i=0;i<n;++i){r[i]=op(i);}
		r.shrk();
		return r;
	}
	template<_mkseq_func F>static poly_t make_i(idt n,F&&op,u32 fx=R2){
		poly_t r{uinit0{n}};
		for(idt i=0;i<n;++i){r[i]=op(i);}
		r.shrk(),r*=fx;
		return r;
	}
	static poly_t make_n(std::initializer_list<u32> lis,u32 fx=R2){
		poly_t r{lis};
		r*=fx;
		return r;
	}
	template<class _it>static poly_t make_s(_it bg,_it ed){
		poly_t r{uinit0(ed-bg)};
		std::copy(bg,ed,r.p),r.shrk();
		return r;
	}
};
using poly=poly_t<aligned_alcor<u32,32> >;
}
#include <sys/mman.h>
#include <sys/stat.h>
namespace QIO_base{
    constexpr int O_buffer_default_size = 1 << 18;
	constexpr int O_buffer_default_flush_threshold = 40;
	constexpr u64 E16 = 1e16, E12 = 1e12;
	constexpr u32 E8 = 1e8, E4 = 1e4;
	struct ict{
		int num[10000];
		constexpr ict(){
			int j = 0;
			for(int e0 = (48 << 0); e0 < (58 << 0); e0 += (1 << 0)){
				for(int e1 = (48 << 8); e1 < (58 << 8); e1 += (1 << 8)){
					for(int e2 = (48 << 16); e2 < (58 << 16); e2 += (1 << 16)){
						for(int e3 = (48 << 24); e3 < (58 << 24); e3 += (1 << 24)){
							num[j] = e0 ^ e1 ^ e2 ^ e3, ++j;
						}
					}
				}
			}
		}
	}constexpr ot;
}
namespace QIO_I {
	using namespace QIO_base;
	struct Qinf{
		FILE* f;
		char *bg,*ed,*p;
		struct stat Fl;
		Qinf(FILE *fi) : f(fi){
			int fd = fileno(f);
			fstat(fd, &Fl);
			bg = (char*)mmap(0, Fl.st_size + 1, PROT_READ,MAP_PRIVATE, fd, 0);
			p = bg, ed = bg + Fl.st_size;
		}
		~Qinf(){
			munmap(bg,Fl.st_size + 1);
		}
		void skip_space(){
			while(*p <= ' '){
				++p;
			}
		}
		char get(){
			return *p++;
		}
		char seek()const{
			return *p;
		}
		Qinf& read(char* s, size_t count){
			return memcpy(s, p, count), p += count, *this;
		}
		template<std::unsigned_integral T>Qinf& operator >> (T&x){
			skip_space(),x=0;
			for(; *p > ' '; ++p){
				x = x * 10 + (*p & 0xf);
			}
			return *this;
		}
	}qin(stdin);
}
namespace QIO_O{
	using namespace QIO_base;
	struct Qoutf{
		FILE *f;
		char *bg,*ed,*p;
		char *ed_thre;
		Qoutf(FILE *fo,size_t sz = O_buffer_default_size):
			f(fo),
			bg(new char[sz]),ed(bg+sz),p(bg),
			ed_thre(ed - O_buffer_default_flush_threshold){}
		void flush(){
			fwrite_unlocked(bg,1,p - bg,f),p=bg;
		}
		void chk(){
			if(__builtin_expect(p > ed_thre, 0)){
				flush();
			}
		}
		~Qoutf(){
			flush();
			delete[] bg;
		}
		void put4(u32 x) {
			auto C = (const char*)(ot.num + x);
			if (x > 99u) {
				if (x > 999u){
					memcpy(p, C, 4), p += 4;
				}
				else{
					memcpy(p, C + 1, 3), p += 3;
				}	
			} 
			else {
				if (x > 9u){
					memcpy(p, C + 2, 2), p += 2;
				}
				else{
					*p++ = x ^ 48;
				}
			}
		}
		void put2(u32 x) {
			if (x > 9u){	
				memcpy(p, (const char*)(ot.num + x) + 2, 2), p += 2;
			}
			else{
				*p++ = x ^ 48;
			}
		}
		Qoutf &write(const char *s, size_t count) {
			if (count > 1024 || p + count > ed_thre)
				flush(), fwrite_unlocked(s, 1, count, f);
			else
				memcpy(p, s, count), p += count, chk();

			return *this;
		}
		Qoutf &operator << (char ch) {
			return *p++ = ch, *this;
		}
		Qoutf &operator << (u32 x) {
			if (x >= E8) {
				put2(x / E8), x %= E8;
				memcpy(p, ot.num + x / E4, 4), p += 4;
				memcpy(p, ot.num + x % E4, 4), p += 4;
			} else if (x >= E4) {
				put4(x / E4);
				memcpy(p, ot.num + x % E4, 4), p += 4;
			} else{
				put4(x);
			}
			return chk(), *this;
		}
	}qout(stdout);
}
namespace QIO{
	using QIO_I::Qinf;
	using QIO_I::qin;
	using QIO_O::Qoutf;
	using QIO_O::qout;
}
using namespace QIO;
void solve(){
	using namespace No_Poly;
	idt N,totd=0;
    qin>>N;
    auto prod=[&](auto&&self,idt l,idt r)->poly{
        auto len=r-l;
        if(len==0){return {R};}
        if(len==1){
            idt d;
            qin>>d,totd+=d;
            poly res;
            res.read(d+1,qin);
            return res;
        }
        auto mid=(l+r)>>1;
        auto lhs=self(self,l,mid);
        lhs*=self(self,mid,r);
        return lhs;
    };
    auto res=prod(prod,0,N);
    res.write(totd+1,qout);
}
int main(){
	std::ios::sync_with_stdio(false);
	std::cin.tie(nullptr);
	solve();
	return 0;
}
/*
5
1 1 5
1 1 4
1 1 3
1 1 2
1 1 1
*/