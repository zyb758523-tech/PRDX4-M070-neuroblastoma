# GitHub 公共仓库与 Zenodo DOI 发布步骤

## 一、发布前先确认

1. 四位作者英文姓名和顺序：Yanbo Zeng、Jinghao Yan、Shuixue Li、Ling Zhou。
2. 四位作者单位是否均为 `Xinjiang Uygur Autonomous Region Children's Hospital`。
3. 是否同意代码使用 MIT License。MIT 允许他人使用、修改和再发布代码，但必须保留版权和许可声明。
4. 如作者有 ORCID，可在正式发布前加入 `CITATION.cff` 和 `.zenodo.json`；ORCID 不是生成 DOI 的必需项。
5. 确认仓库不包含账号密码、访问令牌、非公开患者身份信息或无权再分发的数据。

## 二、建立 GitHub 仓库

推荐仓库名称：

`PRDX4-M070-neuroblastoma`

建议简介：

`Audited R code and frozen M070 model for a single-cell lactate-state and externally validated prognostic study in neuroblastoma.`

操作：

1. 登录 GitHub。
2. 点击右上角 `+`，选择 `New repository`。
3. Repository name 填 `PRDX4-M070-neuroblastoma`。
4. Description 填上面的英文简介。
5. 选择 `Public`。
6. 不要勾选自动生成 README、.gitignore 或 License，因为本包已经包含。
7. 点击 `Create repository`。
8. 将本文件夹中的全部内容上传到仓库根目录，而不是把最外层文件夹或 ZIP 作为单个文件上传。
9. 提交说明填写：`Prepare audited reproducibility release v1.0.0`。

## 三、补回 GitHub 地址

仓库建立后，把下列占位符：

`__GITHUB_REPOSITORY_URL__`

统一替换为真实地址，例如：

`https://github.com/USERNAME/PRDX4-M070-neuroblastoma`

需要替换的文件：

- `CITATION.cff`
- `CODE_AVAILABILITY.md`

提交说明：`Add public repository URL`。

## 四、连接 Zenodo

1. 打开 https://zenodo.org/ 并登录。
2. 建议使用 GitHub 或 ORCID 登录。
3. 在 Zenodo 头像菜单中进入 `GitHub`。
4. 点击 `Sync now`。
5. 找到 `PRDX4-M070-neuroblastoma`，打开右侧开关。

## 五、在 GitHub 发布 v1.0.0

1. 回到 GitHub 仓库主页。
2. 点击右侧 `Releases`，再点击 `Draft a new release`。
3. Tag 填 `v1.0.0`。
4. Release title 填 `PRDX4-M070 neuroblastoma analysis v1.0.0`。
5. Release description 可填：

   `Audited manuscript-associated release containing 11 R scripts, the frozen M070 model, stage result archives, selected outputs, session records, input manifest, and cryptographic checksums.`

6. 点击 `Publish release`。
7. Zenodo 会自动接收该 release 并创建软件归档。

## 六、获得 DOI 后最后回填

Zenodo 完成处理后会为 v1.0.0 生成可永久引用的 DOI。论文应引用与实际分析代码完全对应的这一版本 DOI。以后若发布新版本，Zenodo 会为新版本建立单独的持久记录，并将各版本关联起来。

将版本 DOI 的完整链接，例如 `https://doi.org/10.5281/zenodo.xxxxxxx`，替换：

`__ZENODO_VERSION_DOI_URL__`

需要更新：

- `CODE_AVAILABILITY.md`
- 正文的 Availability of data and materials / Code availability
- 投稿系统中要求填写的代码或数据仓库链接

同时在 `CITATION.cff` 中加入：

`doi: 10.5281/zenodo.xxxxxxx`

并将 `repository-code` 保留为 GitHub 地址。

## 七、DOI 已生成后不要做什么

- 不要删除 v1.0.0 release。
- 不要用新结果覆盖 v1.0.0。
- 如需修改代码，发布 v1.0.1 或 v1.1.0，并让 Zenodo 生成新的版本 DOI。
- 投稿正文必须引用与实际分析一致的版本 DOI，而不是尚未发布的草稿地址。
