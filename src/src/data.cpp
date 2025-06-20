#include "data.hpp"

Data::Data(int n_, int nVar_) {
  n = n_;
  nVar = nVar_;
  data.resize(nVar * n, 0.0);
}

void Data::solInit(int n_, int nvar_) {
  n = n_;
  nVar = nvar_;
  data.resize(nVar * n, 0.0);

  for (int i = 0; i < n; i++) {
    real h = 2.0 / n;
    real xi = h / 2.0 + i * h - 1.0;
    // data[i]=-sin(M_PI*xi);//for burgers equation

    // for sod tube 1D

    real gamma = GAMMA;
    if (xi < 0) {
      (*this)(i, 0) = 1.0;
      (*this)(i, 1) = 0;
      (*this)(i, 2) = 1.0 / (gamma - 1) * 1;
    } else {
      (*this)(i, 0) = 0.125;
      (*this)(i, 1) = 0;
      (*this)(i, 2) = 1.0 / (gamma - 1) * 0.1;
    }
  }
}
void Data::init(int n_, int nvar_) {
  n = n_;
  nVar = nvar_;
  data.resize(nVar * n, 0.0);
}

real &Data::operator()(int i, int ivar) { return data[i * nVar + ivar]; }

real &Data::operator[](int i) {
  // return (this->data)[i];
  return data.at(i);
}

real Data::maxElement(int ivar) {
  if (data.empty())
    return 0;
  real res = data[0];
  for (int i = 1; i < n; i++) {
    res = (res < data[i]) ? data[i] : res;
  }
  return res;
}

void Data::setValue(std::vector<real> value) {
  assert(value.size() == data.size());
  std::copy(value.begin(), value.end(), data.begin());
}

// void Data::operator= (Data& dat)
// {
//     if(this->n==dat.n&&this->nVar==dat.nVar)
//     {
//         for(ind i = 0; i < n*nVar; i++)
//         {
//             (*this)[i]=dat[i];
//         }

//     }
//     else
//     {
//         std::cout<<"incorrect size at class Data=Data\n";
//     }
// }

void Data::operator+=(std::vector<real> arr) {
  assert(arr.size() == n * nVar);
  for (int i = 0; i < n * nVar; i++) {
    (*this)[i] += arr[i];
  }
}

void Data::setZeros() { std::fill(data.begin(), data.end(), 0.0); }

std::vector<real>::iterator Data::begin() { return data.begin(); }
std::vector<real>::iterator Data::random(int i) {
  return data.begin() + i * nVar;
}
std::vector<real>::iterator Data::end() { return data.end(); }

int Data::size() { return data.size(); }

int Data::getNVar() { return nVar; }
int Data::getN() { return n; }

std::vector<real> Data::getIVar(int ivar) {
  std::vector<real> res;
  assert(ivar <= nVar);

  res.reserve(n);
  for (int i = 0; i < n; i++) {
    res.push_back((*this)(i, ivar));
  }
  return res;
}

Data::Data(const Data &origin_) {
  n = origin_.n;
  nVar = origin_.nVar;
  data.resize(n * nVar);
  std::copy(origin_.data.begin(), origin_.data.end(), data.begin());
}

void Data::setvarName(std::vector<std::string> &&varName_) {
  assert(nVar == varName_.size());
  varName.resize(nVar);
  std::copy(varName_.begin(), varName_.end(), varName.begin());
}

real Data::getL2(int ivar) {
  if (ivar >= nVar) {
    std::cout << "Data error: ivar >= nVar";
    return 0;
  }
  real res = 0;
  for (int i = 0; i < n; i++) {
    res += pow((*this)(i, ivar), 2);
  }
  res = std::sqrt(res / n);
  return res;
}
real Data::getLinf(int ivar) {
  if (ivar >= nVar) {
    std::cout << "Data error: ivar >= nVar";
    return 0;
  }
  real res = std::abs((*this)(0, ivar));
  for (int i = 1; i < n; i++) {
    if (res < std::abs((*this)(i, ivar)))
      res = std::abs((*this)(i, ivar));
  }
  return res;
}