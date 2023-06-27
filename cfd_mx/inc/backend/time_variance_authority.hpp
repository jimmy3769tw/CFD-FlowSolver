#pragma once
#include <vector>
// https://en.wikipedia.org/wiki/Time_Variance_Authority for Name of class

class TimeVarianceAuthority {
 public:
  auto GetTime() { return time_; }
  auto SetTerminalTime(double terminal_time) {
    terminal_time_ = terminal_time;
    loop_max_ = terminal_time_ / dt_;
  }

  TimeVarianceAuthority& SetDt(double t) {
    dt_ = t;
    IniTimeRecoding();
    return *this;
  }

  TimeVarianceAuthority& SetWritingDt(double writing_dt) {
    writing_dt_ = writing_dt;
    return *this;
  }

  double GetPercentageOfWriting() { return (time_ / writing_dt_); }
  auto GetFileNo() { return cnt_file_; }
  bool IsWritingTime() {
    auto time_file_next = writing_time_ + writing_dt_;
    if (time_ >= time_file_next) {
      cnt_file_++;
      writing_time_ += writing_dt_;
      return true;
    } else {
      return false;
    }
  }
  void AddLoop() {
    time_ += dt_;
    loop_++;
    dt_ = dt_next_;
  }

  bool IsLoopFinish() {
    return time_ > terminal_time_+dt_;
  }

  TimeVarianceAuthority& SetTime(double t) {
    time_ = t;
    return *this;
  }

  auto GetDt() { return dt_;}
  auto GetLoop(){return loop_;}

 private:
  double dt_;                  // curr dt
  std::vector<double> dt_pre;  // previous dt
  double dt_init_, dt_next_;
  double terminal_time_;
  double time_ = 0.0;
  double writing_dt_;
  double writing_time_;
  size_t loop_ = 1;
  size_t cnt_file_ = 0;
  size_t loop_max_;

  void IniTimeRecoding() {
    dt_init_ = dt_;
    dt_next_ = dt_;
    dt_pre.resize(4, dt_);
  }

  void RecordDt() {
    for (size_t j = dt_pre.size() - 1; j > 0; --j) {
      dt_pre[j] = dt_pre[j - 1];
    }
    dt_pre[0] = dt_;
  }
};