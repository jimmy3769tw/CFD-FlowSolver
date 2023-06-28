
namespace solver{
    using namespace std;

    template<typename matrixT>
    class cg
    {
        public:

        cg(matrixT & lhs_mat){ init(lhs_mat);}

        virtual ~cg(){ }

        // ! init
        inline bool init(const matrixT & lhs_mat)
        {
            lhs_mat_ = lhs_mat;
            length_ = lhs_mat.row();
            p_.resize(length_);
            r_.resize(length_);
            r2_.resize(length_);
            v_.resize(length_);
            ss_.resize(length_);
            t_.resize(length_);

            return true;
        }

        void SetTolerance(double tolerance) { zeta_ = tolerance;}

        std::pair<int, double> solve(const vector<double> & rhs, vector<double> & x); 

        private:
            //cg variables
            int length_;
            double zeta_{1e-5};
            int iters_max_{3000};
            double  alpha_, beta_, norm0_, 
                    norm_, sum_, scal_, 
                    norm1_, norm2_,
                    omega_, rho1_, rho2_;

            vector<double> p_;
            vector<double> r_;
            vector<double> r2_;
            vector<double> v_;
            vector<double> ss_;
            vector<double> t_;

            matrixT lhs_mat_;

        inline double inner_product(const vector<double> & a, const vector<double> & b)
        {
            double r{0.0};

            #pragma omp simd reduction(+:r)
            for(int i=0; i<length_; ++i)
            {  r += a[i] * b[i]; }

            return r;
        }
        inline double minus(const vector<double> & a, const vector<double> & b, vector<dobule> r){

            for(int i=0; i<length_; ++i)
            {  r = a[i] - b[i]; }

        }

        inline void copy(const vector<double> & a, vector<double> & r){
            #pragma omp parallel for
            for(int i=0; i< a.size() ; ++i){
                r[i] = a[i];
            }
        }

    };


    // ! main
    template<typename matrixT>
    inline std::pair<int, double> 
    cg<matrixT>::solve(
        const vector<double> & rhs, vector<double> & x)
    {

        lhs_mat_.multiply(x, p_);

        std::vector<double> r_(m);

        minus(B, p_, r_)

        copy(r_, p_);

        auto nu = inner_product(r_, r_);
        auto norm0 = inner_product(B ,B);

        const auto tolerance = std::pow((zeta_ * norm0), 2);

        int iter= 1; double error;

        for (double beta = 0.0; iter < iters_max_ ;++iter){


            auto q = multiply(p_);

            auto alpha =  nu / std::inner_product(p_.begin(), p_.end(),q.begin() ,0.0);


            for (int i = 0 ; i < m ;++i )
                x[i] += alpha * p_[i];

            for (int i = 0 ; i < m ;++i )
                r_[i] -= alpha * q[i];


            error = std::inner_product(r_.begin(), r_.end(), r_.begin() ,0.0);
            if (error < tolerance){
                break;
            }

            auto mu = std::inner_product(r_.begin(), r_.end(), r_.begin() ,0.0);

            beta = mu / nu;

            for (int i = 0 ; i < m ;++i ) {
                p_[i] = r_[i] + beta *p_[i];
            }

            nu = mu;
        }
        return std::make_pair( iter, std::sqrt(iter) / norm0);
    }
}  // namespace cg



template<typename T>
inline std::pair<int, double> ELL_matrix<T>::npc_cg(
    std::vector <double> &B,
    std::vector <double> &x
){

}
