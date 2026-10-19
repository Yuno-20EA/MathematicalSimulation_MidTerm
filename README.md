# Steepest Descent vs Newton + Backtracking Line Search trên hàm Rosenbrock

Cài đặt bằng **MATLAB / GNU Octave** (không cần toolbox) hai phương pháp tối ưu không ràng buộc —
**phương pháp hướng dốc nhất (steepest descent)** và **phương pháp Newton** — kết hợp
**tìm kiếm đường đi lùi (backtracking line search, Thuật toán 3.1)** để cực tiểu hóa hàm Rosenbrock (3.27).
Chiều dài bước ban đầu $s_0 = 1$ và chiều dài bước được in ra ở **mỗi vòng lặp**.
Thử với hai điểm xuất phát: $x_0=(1.2,\,1.2)$ và $x_0=(-1.2,\,1)$.

## Mục lục
1. [Cấu trúc repo & cách chạy](#1-cấu-trúc-repo--cách-chạy)
2. [Công thức, giải thích và lý do sử dụng](#2-công-thức-giải-thích-và-lý-do-sử-dụng)
3. [Giải thích thuật toán](#3-giải-thích-thuật-toán)
4. [Đầu ra mong muốn](#4-đầu-ra-mong-muốn)
5. [Mục đích & ứng dụng](#5-mục-đích--ứng-dụng)
6. [Kết quả](#6-kết-quả) (mỗi mục: so với mong muốn → hình → đủ yêu cầu chưa → cải tiến)
7. [Tổng kết & hướng phát triển](#7-tổng-kết--hướng-phát-triển)
8. [Tài liệu tham khảo](#8-tài-liệu-tham-khảo)

---

## 1. Cấu trúc repo & cách chạy

```
.
├── run_demo.m                 # điểm vào: chạy toàn bộ, in bảng, vẽ 3 hình, lưu vào results/
├── src/
│   ├── rosenbrock.m           # f, gradient, Hessian
│   ├── backtracking.m         # Thuật toán 3.1 (điều kiện Armijo)
│   ├── line_search_descent.m  # vòng lặp chung cho steepest descent & Newton
│   ├── default_options.m      # mọi tham số (s0, rho, c, tol, điểm xuất phát, ...)
│   ├── print_tables.m         # bảng Iteration | f(x) | Step Length
│   ├── plot_surface3d.m       # Hình 1: mặt 3D + 4 quỹ đạo + chú thích
│   ├── plot_zoom.m            # Hình 2: thu phóng từng quỹ đạo
│   ├── plot_convergence.m     # Hình 3: f(x_k), ||grad||, step length theo vòng lặp
│   ├── floor_contours.m, sample_indices.m   # hàm phụ trợ vẽ
├── tests/run_tests.m          # 8 kiểm thử tự động (không mở hình)
├── results/
│   ├── console_output.txt     # log console do run_demo ghi lại
│   └── images/                # fig1_3d_surface.png, fig2_zoom.png, fig3_convergence.png, ...
└── .github/workflows/tests.yml  # CI: chạy run_tests bằng Octave mỗi lần push
```

```matlab
>> run_demo      % chạy toàn bộ (MATLAB hoặc Octave), mất vài giây
>> addpath tests; run_tests    % kiểm thử
```

Mọi tham số chỉnh trong `src/default_options.m`: `s0 = 1`, `rho = 0.5`, `c = 1e-4`, `tol = 1e-6`, `maxit = 50000`,
`printMax` (số vòng lặp đầu được in), `stepDecimals`, ...

> **Lưu ý về ảnh minh họa:** các ảnh trong `results/images/` được tạo bằng Octave. Hãy chụp lại từ MATLAB của bạn
> (đặc biệt là ảnh ghép Figure + Command Window) và ghi đè lên cùng tên file để README hiển thị ảnh của riêng bạn.

---

## 2. Công thức, giải thích và lý do sử dụng

### 2.1 Hàm Rosenbrock (3.27)

$$f(x)=100\,(x_2-x_1^2)^2+(1-x_1)^2,\qquad x^*=(1,1),\quad f(x^*)=0$$

**Vì sao dùng:** đây là bài kiểm tra kinh điển cho thuật toán tối ưu. Hàm có một *thung lũng cong hình quả chuối*
($x_2=x_1^2$) rất hẹp và rất phẳng dọc theo đáy. Tìm ra thung lũng thì dễ, đi dọc đáy để tới $x^*$ mới khó.
Hàm là tổng hai bình phương nên $f\ge 0$ và $x^*$ là cực tiểu toàn cục duy nhất.

### 2.2 Gradient và Hessian

$$\nabla f(x)=\begin{pmatrix}-400x_1(x_2-x_1^2)-2(1-x_1)\\[2pt] 200(x_2-x_1^2)\end{pmatrix},\qquad
\nabla^2 f(x)=\begin{pmatrix}1200x_1^2-400x_2+2 & -400x_1\\ -400x_1 & 200\end{pmatrix}$$

**Vì sao dùng:** gradient cho thông tin bậc 1 (hướng tăng nhanh nhất), Hessian cho thông tin bậc 2 (độ cong).
Điều kiện dừng là $\|\nabla f(x_k)\|\le 10^{-6}$ vì tại cực tiểu trơn thì $\nabla f(x^*)=0$.

Tại $x^*$: $\nabla^2 f(x^*)=\begin{pmatrix}802&-400\\-400&200\end{pmatrix}$ có trị riêng
$\lambda_{\min}\approx 0.3994,\ \lambda_{\max}\approx 1001.6$, nên **số điều kiện**
$\kappa=\lambda_{\max}/\lambda_{\min}\approx 2508$. Đây là nguyên nhân gốc khiến steepest descent chậm (mục 2.5).

### 2.3 Điều kiện Armijo và backtracking (Thuật toán 3.1)

Cho hướng giảm $p_k$ (tức $\nabla f(x_k)^\top p_k<0$), chấp nhận chiều dài bước $s$ khi

$$f(x_k+s\,p_k)\;\le\;f(x_k)+c\,s\,\nabla f(x_k)^\top p_k,\qquad c\in(0,1)$$

**Ý nghĩa:** mức giảm thực tế phải đạt ít nhất tỉ lệ $c$ của mức giảm mà xấp xỉ tuyến tính dự đoán.
Nếu không thỏa, nhân $s\leftarrow\rho s$ ($\rho\in(0,1)$) và thử lại.
**Vì sao dùng:** nếu chỉ cần "$f$ giảm" thì có thể giảm quá ít (hoặc dao động quanh cực tiểu) và không hội tụ;
Armijo đảm bảo giảm *đủ* ở mỗi bước. Dùng $c=10^{-4}$ (rất nhỏ, chỉ loại các bước tệ) và $\rho=0.5$ (cân bằng giữa số lần thử và độ chính xác của bước).
Bắt đầu từ $s_0=1$ vì với Newton đó là bước "tự nhiên" của mô hình bậc hai.

### 2.4 Hai hướng tìm kiếm

| Phương pháp | Hướng | Xuất phát từ | Hướng này luôn giảm? |
|---|---|---|---|
| Steepest descent | $p_k=-\nabla f(x_k)$ | cực tiểu mô hình **tuyến tính** | Có (khi $\nabla f\ne0$) |
| Newton | $\nabla^2 f(x_k)\,p_k=-\nabla f(x_k)$ | cực tiểu mô hình **bậc hai** $f+g^\top p+\tfrac12 p^\top H p$ | Chỉ khi $H\succ0$ |

Với Newton, $g^\top p=-g^\top H^{-1}g<0$ chỉ khi $H$ xác định dương. Rosenbrock có vùng Hessian không xác định dương,
nên code có **lưới an toàn**: nếu $g^\top p\ge0$ (hoặc giải hệ thất bại) thì thay $p_k$ bằng $-\nabla f(x_k)$.
Hệ $Hp=-g$ được giải bằng toán tử `\`, không tính `inv(H)`.

### 2.5 Tốc độ hội tụ (lý do kết quả khác nhau cỡ nghìn lần)

- **Steepest descent:** hội tụ tuyến tính, với bài toán toàn phương và line search chính xác thì
  $f_{k+1}-f^*\le\big(\tfrac{\kappa-1}{\kappa+1}\big)^2(f_k-f^*)$. Với $\kappa\approx2508$ hệ số này $\approx0.9984$,
  tức muốn giảm $f$ thêm $10^{10}$ lần cần cỡ $\ln(10^{10})/0.0016\approx1.4\times10^{4}$ vòng. Con số này (ước lượng cục bộ) cùng bậc độ lớn với thực nghiệm (~13 000).
- **Newton:** gần $x^*$ với $s=1$ hội tụ **bậc hai**: $\|x_{k+1}-x^*\|\le C\|x_k-x^*\|^2$, số chữ số đúng xấp xỉ nhân đôi mỗi vòng.
  Điều kiện Armijo với $c<\tfrac12$ chấp nhận $s=1$ khi đủ gần nghiệm, nên backtracking không phá tính chất này.

---

## 3. Giải thích thuật toán

**Thuật toán 3.1 – Backtracking line search**
```
Cho x, hướng giảm p, s0 > 0, rho ∈ (0,1), c ∈ (0,1)
s ← s0
trong khi f(x + s p) > f(x) + c · s · ∇f(x)ᵀp :
    s ← rho · s
trả về s
```

**Vòng lặp chung (steepest descent hoặc Newton)**
```
x ← x0
trong khi ||∇f(x)|| > tol và k < maxit :
    p ← −∇f(x)                         (steepest descent)
    p ← −H⁻¹∇f(x), nếu ∇fᵀp ≥ 0 thì p ← −∇f(x)   (Newton, có lưới an toàn)
    s ← backtracking(x, p)             (Thuật toán 3.1, s0 = 1)
    x ← x + s p ;  k ← k + 1 ;  in (k, f(x), s)
```

Hai phương pháp dùng **chung một vòng lặp** (`line_search_descent.m`) và chỉ khác hướng $p_k$, nên mọi khác biệt về kết quả quy được về hướng tìm kiếm.

---

## 4. Đầu ra mong muốn

**Định dạng console** (mỗi điểm xuất phát, mỗi phương pháp một bảng):
```
--- STARTING FROM POINT [1.2, 1.2] ---

Steepest Descent (Backtracking Line Search)
Iteration      f(x)   |   Step Length
        1 |    0.4310 | 0.000977
        ...
Newton's Method (Backtracking Line Search)
...
Newton converged in fewer iterations.
```

**Hình:** (1) mặt 3D + 4 quỹ đạo + chú thích *Start points / Backtracking occurring here / Ideal – [1,1]*;
(2) thu phóng từng quỹ đạo; (3) hội tụ theo vòng lặp.

**Tiêu chí đạt (dùng để chấm ở mục 6):**

| Mã | Tiêu chí |
|---|---|
| E1 | Cả hai phương pháp, cả hai điểm xuất phát, hội tụ tới $(1,1)$ với $\|\nabla f\|\le10^{-6}$ |
| E2 | $f(x_k)$ không tăng; mọi $s_k\in(0,1]$ và có dạng $\rho^m$ |
| E3 | Newton cần **ít vòng lặp hơn rất nhiều** so với steepest descent |
| E4 | Newton chủ yếu nhận bước đầy đủ $s=1$; chỉ co bước ở vài vòng đầu |
| E5 | Steepest descent bị co bước ở (hầu như) mọi vòng ($s\ll1$) |
| E6 | Đuôi hội tụ Newton có dạng bậc hai (số chữ số đúng tăng nhanh) |
| E7 | Điểm khó $(-1.2,1)$ khiến Newton tốn thêm vòng lặp / co bước nhiều hơn điểm dễ |
| E8 | Có chiều dài bước in ở **mỗi** vòng lặp, đúng bố cục yêu cầu |

> **Lưu ý:** ảnh mẫu của bài yêu cầu chỉ để tham khảo bố cục. Một số con số trong ảnh mẫu (ví dụ steepest descent có $s=1$ ở vòng 1,
> hay $f(x)<0$) không thể xảy ra với hàm Rosenbrock thật nên **không** dùng làm đáp án; đáp án là số do code sinh ra.
> Ảnh mẫu ghi điểm xuất phát thứ hai là $[-1.2,\,1.2]$, còn đề bài là $(-1.2,\,1)$; repo theo đề bài.

---

## 5. Mục đích & ứng dụng

**Mục đích học thuật:** thấy tận mắt (bằng số liệu và quỹ đạo) cái giá của việc chỉ dùng thông tin bậc 1 so với bậc 2,
vai trò của line search trong việc bảo đảm hội tụ, và vì sao số điều kiện quyết định tốc độ.

**Ứng dụng của cùng bộ khung "hướng giảm + line search":**
- **Huấn luyện mô hình học máy:** gradient descent là nền của mọi trình tối ưu; các biến thể bậc 2 xấp xỉ (L-BFGS, Gauss–Newton) dùng cho hồi quy logistic và các bài toán cỡ nhỏ–vừa.
- **Thị giác máy tính:** tinh chỉnh tham số camera, ước lượng tư thế (PnP), bundle adjustment đều là bài toán bình phương tối thiểu phi tuyến, thường giải bằng Gauss–Newton hoặc Levenberg–Marquardt, là họ hàng trực tiếp của phương pháp Newton ở đây.
- **Edge AI / hệ nhúng:** hiệu chuẩn cảm biến và tinh chỉnh tham số trên thiết bị với số biến nhỏ, nơi Hessian nhỏ có thể tính và giải được rẻ; bài toán đánh đổi *số vòng lặp* và *chi phí mỗi vòng* chính là bài toán chọn thuật toán cho thiết bị hạn chế bộ nhớ và năng lượng.
- **Vật lý / khoa học dữ liệu:** khớp mô hình vào dữ liệu thực nghiệm, tìm cực tiểu năng lượng.

---

## 6. Kết quả

Tất cả số liệu dưới đây lấy từ lần chạy với `s0=1, rho=0.5, c=1e-4, tol=1e-6` (xem `results/console_output.txt`).

### 6.1 Console – điểm dễ $x_0=(1.2,\,1.2)$

| | Steepest descent | Newton |
|---|---|---|
| Số vòng lặp | 13 011 | **8** |
| Nghiệm cuối | (1.000001, 1.000002) | (1.000000, 1.000000) |
| $f$ cuối / $\|\nabla f\|$ cuối | $6.2\times10^{-13}$ / $9.99\times10^{-7}$ | $1.1\times10^{-25}$ / $1.4\times10^{-11}$ |
| Chiều dài bước | luôn $<1$ (nhỏ nhất $9.77\times10^{-4}$) | `[1 0.5 1 1 1 1 1 1]` |
| Số lần tính $f$ trong backtracking | 129 262 | 9 |
| Lần thay hướng Newton bằng $-\nabla f$ | – | 0 |

Vòng 1 của steepest descent: $f$ giảm từ 5.8 xuống 0.431 với $s\approx9.8\times10^{-4}$; Newton: $f$ giảm xuống 0.0384 với $s=1$.

**Ảnh minh họa:** ![console 1](results/images/console_start1.png) *(chụp màn hình Command Window của bạn, để trống nếu chưa có)*

**So với mong muốn:** E1 ✔, E2 ✔, E3 ✔ (khoảng 1 600 lần ít vòng hơn), E4 ✔ (7/8 bước là $s=1$), E5 ✔ (13 011/13 011 bước $s<1$),
E6 ✔ ($f$: $1.9\times10^{-5}\to3.4\times10^{-8}\to3.2\times10^{-14}\to1.1\times10^{-25}$, số mũ gần nhân đôi), E8 ✔. **Chuẩn xác.**

**Chất lượng so với yêu cầu ban đầu:** đủ. Cả hai hội tụ đúng nghiệm giải tích $(1,1)$; bảng đúng bố cục yêu cầu.

**Có thể cải tiến:** (i) với steepest descent, thay Armijo thuần bằng điều kiện Wolfe hoặc bước Barzilai–Borwein để tránh bước nhỏ cỡ $10^{-3}$;
(ii) thêm tiêu chí dừng theo $f$ và theo $\|x_{k+1}-x_k\|$; (iii) in thêm $\|\nabla f\|$ mỗi vòng.

### 6.2 Console – điểm khó $x_0=(-1.2,\,1)$

| | Steepest descent | Newton |
|---|---|---|
| Số vòng lặp | 13 756 | **21** |
| Nghiệm cuối | (0.999999, 0.999998) | (1.000000, 1.000000) |
| $f$ cuối / $\|\nabla f\|$ cuối | $6.1\times10^{-13}$ / $9.83\times10^{-7}$ | $3.7\times10^{-21}$ / $4.5\times10^{-10}$ |
| Chiều dài bước | luôn $<1$; chỉ 1 bước $>0.01$ | `1, 0.125, 1, 1, 1, 0.25, 1, 1, 1, 0.5, 1, 1, 1, 0.5, 1 ...` (4 bước bị co, 17 bước $s=1$) |
| Số lần tính $f$ trong backtracking | 136 799 | 28 |
| Lần thay hướng Newton bằng $-\nabla f$ | – | 0 |

**Ảnh minh họa:** ![console 2](results/images/console_start2.png) *(để trống để bạn tự chụp)*

**So với mong muốn:** E1 ✔, E2 ✔, E3 ✔ (khoảng 650 lần), E4 ✔ (17/21 bước $s=1$), E5 ✔, E6 ✔, **E7 ✔**
(Newton tốn 21 vòng và 4 lần co bước so với 8 vòng và 1 lần co bước ở điểm dễ). **Chuẩn xác.**

**Chất lượng:** đủ. Điểm khó làm Newton chậm đi nhưng vẫn hội tụ nhờ backtracking, đúng ý nghĩa "bảo vệ toàn cục" của line search.
Lưới an toàn không được kích hoạt lần nào: Hessian dương ở mọi điểm lặp trong lần chạy này (không có nghĩa nó vô dụng với điểm xuất phát khác).

**Có thể cải tiến:** (i) thử thêm các điểm xuất phát ngẫu nhiên và báo cáo thống kê (số vòng, tỉ lệ hội tụ), vì kết luận từ hai điểm là yếu về mặt thống kê;
(ii) dùng Newton cải tiến (cộng $\tau I$ vào Hessian khi không xác định dương) thay vì quay về $-\nabla f$; (iii) khảo sát độ nhạy theo $\rho, c$.

### 6.3 Hình 1 – mặt 3D và 4 quỹ đạo

![3D](results/images/fig1_3d_surface.png)

**So với mong muốn:** đủ các thành phần: mặt Rosenbrock có đường đồng mức ở đáy, 4 quỹ đạo (SD xanh dương/xanh lơ, Newton đỏ/xanh lá),
chú thích *Start point*, *Backtracking occurring here*, *Ideal – [1, 1]*. Bề mặt cắt ở $f=200$ và được đặt hơi trong suốt (`FaceAlpha = 0.7`, nếu phần mềm hỗ trợ) để thấy quỹ đạo trong thung lũng.
Mức đạt: **đạt về nội dung**; vị trí chữ ghi chú có thể cần nhích tay tùy phần mềm/độ phân giải (đã kiểm tra trên Octave, chưa kiểm tra bằng ảnh chụp từ MATLAB).

**Chất lượng:** đủ để minh họa vị trí thung lũng và nghiệm, nhưng góc 3D **không phải** nơi tốt nhất để so sánh hai thuật toán (quỹ đạo nhỏ so với bề mặt). Việc so sánh định lượng nằm ở Hình 2 và 3.

**Có thể cải tiến:** thêm `rotate3d`/nhiều góc nhìn; vẽ mũi tên 3D thật bằng `annotation`; nâng độ phân giải lưới; xuất `.fig`/`.gif` xoay.

### 6.4 Hình 2 – thu phóng từng quỹ đạo (A overview – B 15 vòng đầu – C quanh nghiệm)

![zoom](results/images/fig2_zoom.png)

**So với mong muốn:** hàng trên $x_0=(1.2,1.2)$, hàng dưới $x_0=(-1.2,1)$. Cột (B) cho thấy rõ steepest descent nhảy một bước lớn rồi **nhích rất chậm và lệch hai bên đáy thung lũng**,
còn Newton bám sát đáy thung lũng chỉ với vài điểm lặp được đánh số. Cột (C) cho thấy Newton chạm gần $x^*$ chỉ sau vài vòng cuối trong khi steepest descent vẫn đang trượt dọc đáy.
Đúng với yêu cầu "đồ thị phóng to chỉ rõ sự thay đổi của hai thuật toán, nhìn rõ đường đi". **Đạt.**

**Chất lượng:** đủ; nhược điểm là ở cột (C) quỹ đạo hai phương pháp chồng lên nhau gần $x^*$ (vì cùng nằm trong đáy), khó phân biệt bằng mắt.

**Có thể cải tiến:** vẽ cửa sổ (C) với tỉ lệ log theo $\|x-x^*\|$; thêm mũi tên chỉ hướng bước; thêm thanh trượt chọn số vòng lặp hiển thị (`uicontrol`/`slider`); xuất video quỹ đạo từng vòng.

### 6.5 Hình 3 – hội tụ theo vòng lặp

![convergence](results/images/fig3_convergence.png)

**So với mong muốn:** cột trái $f(x_k)$ và cột giữa $\|\nabla f\|$ (log–log): đường Newton rơi dốc xuống trong chục vòng, steepest descent đi gần phẳng rồi giảm chậm trong $\sim10^4$ vòng;
cột phải: chiều dài bước (log): Newton phần lớn bằng 1, steepest descent bám quanh $10^{-3}$. Khớp E3, E4, E5, E6. **Đạt.**

**Chất lượng:** đủ; trục log theo $k$ là cần thiết vì hai phương pháp chênh nhau ba bậc độ lớn về số vòng lặp.

**Có thể cải tiến:** vẽ theo **thời gian chạy** và **số lần tính hàm** thay vì theo vòng lặp (so sánh công bằng vì mỗi vòng Newton đắt hơn: thêm Hessian và giải hệ);
thêm đường tham chiếu lý thuyết $\big(\tfrac{\kappa-1}{\kappa+1}\big)^{2k}$.

### 6.6 Kiểm thử tự động

`run_tests` kiểm tra: gradient và Hessian khớp sai phân hữu hạn; Hessian tại $x^*$ xác định dương với $\kappa\approx2.5\times10^3$;
bước backtracking thỏa Armijo và có dạng $\rho^m$ (và bước gấp đôi vi phạm Armijo); hai phương pháp hội tụ từ hai điểm xuất phát; $f$ giảm đơn điệu và $s\in(0,1]$; Newton ít vòng lặp hơn $\ge100$ lần.
Kết quả: **8/8 đạt** (Octave). **So với mong muốn:** đúng, và là cơ sở để tin các con số ở trên.
**Cải tiến:** thêm test hồi quy so sánh số vòng lặp với giá trị tham chiếu, thêm test cho nhánh lưới an toàn (cần điểm có Hessian không xác định dương), chạy CI cả trên MATLAB.

---

## 7. Tổng kết & hướng phát triển

**Kết luận:** mọi tiêu chí E1–E8 đạt ở cả hai điểm xuất phát. Newton + backtracking hội tụ sau 8 và 21 vòng, steepest descent cần hơn $1.3\times10^4$ vòng.
Sự chênh lệch này khớp với phân tích số điều kiện $\kappa\approx2508$.

**Hạn chế hiện tại:** chỉ 2 điểm xuất phát; một bộ $(c,\rho)$; chưa đo thời gian; chưa có điều kiện Wolfe (nên chưa có đảm bảo cho BFGS).

**Hướng phát triển:**
1. Quasi-Newton (BFGS, L-BFGS) với điều kiện Wolfe: giữ ưu thế của Newton mà không cần Hessian.
2. Trust-region và Levenberg–Marquardt thay cho line search.
3. Thống kê trên nhiều điểm xuất phát; khảo sát độ nhạy $\rho, c, s_0$.
4. Mở rộng lên $n$ chiều (Rosenbrock tổng quát) và so sánh chi phí $O(n)$ với $O(n^3)$ mỗi vòng.
5. **Cài lại bằng C++ (Eigen)** rồi thử trên vi điều khiển, bước nối từ phần toán sang lập trình nhúng.

---

## 8. Tài liệu tham khảo
- J. Nocedal, S. J. Wright, *Numerical Optimization*, 2nd ed., Springer, 2006 (chương 2–3: line search, Armijo, backtracking, Newton).
- H. H. Rosenbrock, "An automatic method for finding the greatest or least value of a function", *The Computer Journal* 3(3), 1960.
- L. Armijo, "Minimization of functions having Lipschitz continuous first partial derivatives", *Pacific J. Math.* 16(1), 1966.

## Giấy phép
MIT (xem `LICENSE`; thay `<TEN-CUA-EM>` bằng tên của bạn).
