# Plasma

Calculation of the plasma scattering spectra

Input plasma parameters see in Plasma_funcs.h. They also listed in Pars.config.

Output table format:

frequency \t amplitude

Usually consider spectrum at frequency band of -5:5 Khz

## Как собирать

Для компиляции проекта требуется компилятор с поддержкой C++14. Проект кроссплатформенный и тестировался только на компиляторе GCC (либо его версиях Clang и MinGW для macOS и Windows соответственно).

### Windows

Нужно установить [CMake](https://cmake.org/download/).

#### Установка компилятора

На [официальном сайте MinGW-w64](https://www.mingw-w64.org/downloads/) есть список различных источников. Я рекомендую качать **"Mingw-builds"** с [GitHub](https://github.com/niXman/mingw-builds-binaries/releases). Качать надо [**x86_64-13.1.0-release-posix-seh-msvcrt-rt_v11-rev1.7z**](https://github.com/niXman/mingw-builds-binaries/releases/download/13.1.0-rt_v11-rev1/x86_64-13.1.0-release-posix-seh-msvcrt-rt_v11-rev1.7z). После скачивания нужно распаковать архив куда угодно и добавить в переменную среды `%PATH%` следующие пути (например, пусть корневой каталог компилятора таков: *C:\\MinGW*):

- C:\MinGW\bin;
- C:\MinGW\libexec\gcc\x86_64-w64-mingw32\13.1.0.

#### Сборка

В каталоге, в котором будет производиться сборка, нужно запустить **CMake**:

```cmd
cmake %path_to_Plasma% -G "MinGW Makefiles" && mingw32-make
```

где `%path_to_Plasma%` — путь до данного репозитория.
